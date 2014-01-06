#include "board.h"
#include "mw.h"

typedef struct fp_vector
{
    float X, Y, Z;
} t_fp_vector_def;

typedef union
{
    float A[3];
    t_fp_vector_def V;
} t_fp_vector;

float        accSmooth[3], ACC_speed[2];
float        accADC[3], gyroADC[3], magADCfloat[3];
int32_t      sonarAlt;
float        BaroAlt, EstAlt, AltHold, vario;            // variometer in cm/s + is up
int16_t      BaroP, BaroI, BaroD;
bool         newbaroalt, GroundAltInitialized;
float        ACCDeltaTimeINS = 0;

// **************
// gyro+acc IMU
// **************
/*
 * Sensor data rate:
 * Baro  - 25 Hz  - 40 ms  | 50 ms
 * Accel - 400 Hz - 2.5 ms | 10 ms
 * Mag   - 30 Hz  - 4.5 ms | 40 ms
 * Gyro  - 760 Hz - 1.3 ms | 10 ms
 */

float   gyroData[3] = { 0, 0, 0 }, angle[2] = { 0, 0 };                    // absolute angle inclination in multiple of 0.1 degree    180 deg = 1800
static  uint8_t SmoothingFactor[3]  = { 0, 0, 0 };
static  bool    GyroSmoothing;

static void getEstimatedAttitude(void);

void imuInit(void)                                                         // Initialize & precalculate some values here
{
    if (cfg.gy_smrll || cfg.gy_smptc || cfg.gy_smyw)
    {
        SmoothingFactor[ROLL]  = cfg.gy_smrll;
        SmoothingFactor[PITCH] = cfg.gy_smptc;
        SmoothingFactor[YAW]   = cfg.gy_smyw;
        GyroSmoothing          = true;
    } else GyroSmoothing = false;

#ifdef MAG
    if (sensors(SENSOR_MAG)) Mag_init();
#endif
}

void computeIMU(void)
{
    static  float LastGyroSmooth[3] = { 0.0f, 0.0f, 0.0f };
    static  int16_t triywavg[4];
    static  uint8_t triywavgpIDX = 0;
    uint8_t axis, i;
    float   flttmp;

    if (MpuSpecial)
    {
        GETMPU6050();
        getEstimatedAttitude();
    }
    else
    {
        gyro.temperature(&telemTemperature1);                    // Read out gyro temperature
        Gyro_getADC();                                           // Also feeds gyroData
        if (sensors(SENSOR_ACC))
        {
            ACC_getADC();
            getEstimatedAttitude();
        }
    }

    if(cfg.mixerConfiguration == MULTITYPE_TRI && cfg.gy_smyw)   // Moving average for yaw in tri mode
    {
        triywavg[triywavgpIDX] = (int16_t)gyroData[YAW]; triywavgpIDX++;
        if (triywavgpIDX == 4) triywavgpIDX = 0;
        flttmp = 0;
        for (i = 0; i < 4; i++) flttmp += triywavg[i];
        gyroData[YAW] = flttmp * 0.25f;
    }

    if (GyroSmoothing)
    {
        for (axis = 0; axis < 3; axis++)
        {
            if (SmoothingFactor[axis] > 1)                       // Circumvent useless action
            {
                flttmp               = (float)SmoothingFactor[axis];
                gyroData[axis]       = ((LastGyroSmooth[axis] * (flttmp - 1.0f)) + gyroData[axis]) / flttmp;
                LastGyroSmooth[axis] = gyroData[axis];
            }
        }
    }
}

// Rotate Estimated vector(s) with small angle approximation, according to the gyro data
void rotateV(struct fp_vector *v, float *delta)
{
    struct    fp_vector v_tmp = *v;
    float     mat[3][3];                                                  // This does a  "proper" matrix rotation using gyro deltas without small-angle approximation
    float     cosx, sinx, cosy, siny, cosz, sinz;
    float     coszcosx, coszcosy, sinzcosx, coszsinx, sinzsinx;
    cosx      = cosf(-delta[PITCH]);
    sinx      = sinf(-delta[PITCH]);
    cosy      = cosf(delta[ROLL]);
    siny      = sinf(delta[ROLL]);
    cosz      = cosf(delta[YAW]);
    sinz      = sinf(delta[YAW]);
    coszcosx  = cosz * cosx;
    coszcosy  = cosz * cosy;
    sinzcosx  = sinz * cosx;
    coszsinx  = sinx * cosz;
    sinzsinx  = sinx * sinz;
    mat[0][0] = coszcosy;
    mat[0][1] = sinz * cosy;
    mat[0][2] = -siny;
    mat[1][0] = (coszsinx * siny) - sinzcosx;
    mat[1][1] = (sinzsinx * siny) + (coszcosx);
    mat[1][2] = cosy * sinx;
    mat[2][0] = (coszcosx * siny) + (sinzsinx);
    mat[2][1] = (sinzcosx * siny) - (coszsinx);
    mat[2][2] = cosy * cosx;
    v->X      = v_tmp.X * mat[0][0] + v_tmp.Y * mat[1][0] + v_tmp.Z * mat[2][0];
    v->Y      = v_tmp.X * mat[0][1] + v_tmp.Y * mat[1][1] + v_tmp.Z * mat[2][1];
    v->Z      = v_tmp.X * mat[0][2] + v_tmp.Y * mat[1][2] + v_tmp.Z * mat[2][2];
}

// *Somehow* modified by me..
// Besides mwii credit must go to Sebbi and BRM! Hopefully they condone mentioning them above my trash.
// Sebbi for his rotation of the acc vector and BRM for his normalization ideas.
static void getEstimatedAttitude(void)
{
    static t_fp_vector EstG, EstM;
    static float    accLPFINS[3];
    static float    INV_GYR_CMPF_FACTOR, INV_GYR_CMPFM_FACTOR, ACC_INS_RC, ACC_RC;
    static uint32_t previousT;
    static bool     init = false;
    float           scale, deltaGyroAngle[3], ACCinsFac, ACCFac, rollRAD, pitchRAD;
    float           cr, sr, cp, sp, cy, sy, spcy, spsy, acc_south, acc_west, acc_up;
    float           tmp[4], AccMag = 0;
    uint8_t         axis;
    uint32_t        currentT = micros();

    tmp[0]          = (float)(currentT - previousT);
    scale           = tmp[0] * GyroScale;
    ACCDeltaTimeINS = tmp[0] * 0.000001f;
    previousT       = currentT;

    if(!init)
    {
        init = true;
        INV_GYR_CMPF_FACTOR  = 1.0f / (float)(cfg.gy_cmpf  + 1);               // Default 400
        INV_GYR_CMPFM_FACTOR = 1.0f / (float)(cfg.gy_cmpfm + 1);               // Default 200
        ACC_INS_RC           = 0.5f / (M_PI * cfg.acc_ilpfhz);                 // Default 5,895 Hz
        ACC_RC               = 0.5f / (M_PI * cfg.acc_lpfhz);                  // Default 0,536 Hz
        for (axis = 0; axis < 3; axis++)                                       // Preset some values to reduce runup time
        {
            accLPFINS[axis] = accADC[axis];
            accSmooth[axis] = accADC[axis];
            if (axis == YAW) EstG.A[axis] = acc_1G;
            else EstG.A[axis] = 0.0f;
        }
    }

    ACCinsFac = ACCDeltaTimeINS / (ACC_INS_RC + ACCDeltaTimeINS);             // Adjust LPF to cycle time / do Hz cut off
    ACCFac    = ACCDeltaTimeINS / (ACC_RC     + ACCDeltaTimeINS);
    for (axis = 0; axis < 3; axis++)
    {
        deltaGyroAngle[axis] = gyroADC[axis]   * scale;
        accLPFINS[axis]     += ACCinsFac       * (accADC[axis] - accLPFINS[axis]);
        accSmooth[axis]     += ACCFac          * (accADC[axis] - accSmooth[axis]);
        AccMag              += accSmooth[axis] * accSmooth[axis];
    }
    AccMag = (AccMag * 100) / SQacc_1G;
    rotateV(&EstG.V, deltaGyroAngle);
    // Apply complimentary filter (Gyro drift correction)
    // If accel magnitude >1.15G or <0.85G and ACC vector outside of the limit range => we neutralize the effect of accelerometers in the angle estimation.
    // To do that, we just skip filter, as EstV already rotated by Gyro
    if (72 < AccMag && AccMag < 133)
    {
        for (axis = 0; axis < 3; axis++) EstG.A[axis] = (EstG.A[axis] * (float)cfg.gy_cmpf + accSmooth[axis]) * INV_GYR_CMPF_FACTOR;
    }

    if (EstG.A[YAW] > ACCZ_25deg) f.SMALL_ANGLES_25 = 1;
    else f.SMALL_ANGLES_25 = 0;

    for (axis = 0; axis < 3; axis++) tmp[axis] = EstG.A[axis];
    tmp[3]       =  sqrtf(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
    if (!tmp[3]) tmp[3] = 1;
    else tmp[3]  = 1.0f / tmp[3];
    for (axis = 0; axis < 3; axis++) tmp[axis] *= tmp[3];
    rollRAD      =  atan2f(tmp[0], tmp[2]);
    pitchRAD     =  asinf(-tmp[1]);                                           // Has to have the "wrong sign" relative to angle[PITCH]
    angle[ROLL]  =  rollRAD  * RADtoDEG10;
    angle[PITCH] = -pitchRAD * RADtoDEG10;
    cr           =  cosf(rollRAD);
    sr           =  sinf(rollRAD);
    cp           =  cosf(pitchRAD);
    sp           =  sinf(pitchRAD);
    TiltValue    =  EstG.V.Z * INVacc_1G;                                     // / acc_1G;
    if (sensors(SENSOR_MAG) && cfg.mag_calibrated)                            // Every mag normalization attempt worsened the PH
    {                                                                         // Tried: Normalizing EstM, Normalizing Xh,Yh (so tmp[0], tmp[1])
        rotateV(&EstM.V, deltaGyroAngle);                                     // Rotate According G vector
        for (axis = 0; axis < 3; axis++) EstM.A[axis] = (EstM.A[axis] * (float)cfg.gy_cmpfm + magADCfloat[axis]) * INV_GYR_CMPFM_FACTOR;
        tmp[0]  = EstM.A[1] * cp + EstM.A[0] * sr * sp + EstM.A[2] * cr * sp;
        tmp[1]  = EstM.A[0] * cr - EstM.A[2] * sr;
        heading = wrap_180(atan2f(-tmp[1], tmp[0]) * RADtoDEG + magneticDeclination); // Get rad to Degree and add declination (without *10 shit) // Wrap to -180 0 +180 Degree        
    } else heading = 0;                                                       // if no mag or not calibrated do bodyframe below
    tmp[0]    = heading * RADX;                                               // Do GPS INS rotate ACC X/Y to earthframe no centrifugal comp. yet
    cy        = cosf(tmp[0]);
    sy        = sinf(tmp[0]);
    cos_yaw_x = cy;                                                           // Store for general use
    sin_yaw_y = sy;                                                           // Store for general use
    spcy      = sp * cy;
    spsy      = sp * sy;
    for (axis = 0; axis < 3; axis++) tmp[axis] = accLPFINS[axis] * INVacc_1G; // / acc_1G Reference to Gravity before rotation
    acc_up    = ((-sp) * tmp[1] + sr * cp * tmp[0] + cp * cr * tmp[2]) - 1;   // -1G That works good for althold
    for (axis = 0; axis < 3; axis++) tmp[axis] = accLPFINS[axis];
    tmp[3]    = sqrtf(EstG.V.X * tmp[0] + EstG.V.Y  * tmp[1] + EstG.V.Z * tmp[2]); // Normalize ACCvector so the gps ins works
    if (!tmp[3]) tmp[3] = 1;                                                  // In that case all tmp must be zero so div by 1 is ok
    else tmp[3] = 1.0f / tmp[3];
    for (axis = 0; axis < 3; axis++) tmp[axis] *= tmp[3];
    acc_south = (cp * cy) * tmp[1] + (sr * spcy - cr * sy) * tmp[0] + ( sr * sy + cr * spcy) * tmp[2];
    acc_west  = (cp * sy) * tmp[1] + (cr * cy + sr * spsy) * tmp[0] + (-sr * cy + cr * spsy) * tmp[2];
    tmp[3]    = 980.665f  * ACCDeltaTimeINS;                                  // vel factor for normalized output tmp3      = (9.80665f * (float)ACCDeltaTime) / 10000.0f;
    tmp[2]    = constrain(TiltValue, 0.5f, 1.0f) * tmp[3];                    // Empirical reduction of hightdrop in forward flight
    if(GroundAltInitialized) vario += acc_up * tmp[2];                        // Positive when moving Up. Just do Vario when Baro completely initialized.
    ACC_speed[LAT] -= acc_south * tmp[3];                                     // Positive when moving North cm/sec when no MAG this is speed to the front
    ACC_speed[LON] -= acc_west  * tmp[3];                                     // Positive when moving East cm/sec when no MAG this is speed to the right
}

#ifdef BARO
///////////////////////////////////////////////
//Crashpilot1000 Mod getEstimatedAltitude ACC//
///////////////////////////////////////////////
#define VarioTabsize 8

void getEstimatedAltitude(void)
{
    static uint8_t  Vidx, IniStep = 0, IniCnt = 0;
    static uint32_t LastBarotime = 0;
    static float    VarioTab[VarioTabsize], AvgHz, LastEstAltBaro, SNRcorrect, SNRavg;
    float           fltemp, EstAltBaro;
    uint32_t        TimeTemp;
    uint8_t         i;

    if (!GroundAltInitialized)
    {
        if (newbaroalt)
        {
            TimeTemp     = micros();
            fltemp       = (float)(TimeTemp - LastBarotime);     // fltemp has deltaT in us
            LastBarotime = TimeTemp;
            switch(IniStep)                                      // Casemachine here for further extension
            {
            case 0:
                IniCnt++;
                if(IniCnt == 50)                                 // Waste 50 Cycles to let things (buffers) settle then ini some vars and proceed
                {
                    for (i = 0; i < VarioTabsize; i++) VarioTab[i] = 0;
                    AvgHz  = EstAlt = GroundAlt = vario = SNRavg = 0;
                    IniCnt = SonarStatus = 0;
                    IniStep++;
                }
                break;
            case 1:
                GroundAlt += BaroAlt;
                AvgHz     += fltemp;
                IniCnt++;
                if (IniCnt == 50)                               // Gather 50 values
                {
                    GroundAlt *= 0.02f;
                    AvgHz      = 50000000.0f / AvgHz;           // Calculate Average Hz here since we skip Baro temp readout every 2nd read
                    GroundAltInitialized = true;  
                }
                break;
            }
        }
    }
    else
    {
        if (SonarStatus) fltemp = sonarAlt;
        switch(SonarStatus)
        {
        case 0:
            SNRavg  = 0;
            IniStep = 0;
            break;
        case 1:
            if (!IniStep)
            {
                IniStep = 1;
                SNRavg  = fltemp;
            }
            else SNRavg += 0.2f * (fltemp - SNRavg);             // Adjust Average during accepttimer (ca. 550ms so ca. 20 cycles)
            SNRcorrect = EstAlt + GroundAlt - SNRavg;            // Calculate baro/sonar displacement on 1st contact
            break;
        case 2:
            if (newbaroalt) BaroAlt = (SNRcorrect + fltemp) * cfg.snr_cf + BaroAlt * (1 - cfg.snr_cf); // Set weight / make transition smoother
            break;
        }

        EstAlt += vario * ACCDeltaTimeINS;
        if (newbaroalt)
        {
            EstAltBaro     = BaroAlt - GroundAlt;
            VarioTab[Vidx] = constrain(EstAltBaro - LastEstAltBaro, -127.0f, 127.0f) * AvgHz; Vidx++;
            if (Vidx == VarioTabsize) Vidx = 0;
            LastEstAltBaro = EstAltBaro;
            fltemp = 0;
            for (i = 0; i < VarioTabsize; i++) fltemp += VarioTab[i];
            fltemp /= (float)VarioTabsize;                       // fltemp = BaroClimbRate in cm/sec // + is up // 27ms * 37 = 999ms
            vario   = vario  * cfg.accz_vcf + fltemp     * (1.0f - cfg.accz_vcf);
            EstAlt  = EstAlt * cfg.accz_acf + EstAltBaro * (1.0f - cfg.accz_acf);
            if (cfg.bar_dbg)
            {
                debug[0] = EstAltBaro * 10;
                debug[1] = EstAlt * 10;
                debug[2] = fltemp;
                debug[3] = vario;
            }
        }
    }
}

void getAltitudePID(void)                                        // I put this out of getEstimatedAltitude seems logical
{
    float ThrAngle;
    ThrAngle = constrain(TiltValue * 100.0f, 0, 100.0f);
    BaroP  = BaroI = BaroD = 0;                                  // Reset the Pid, create something new, or not....
    if (ThrAngle < 40 || TiltValue < 0) return;                  // Don't do BaroPID if copter too tilted
    BaroP  = (int16_t)((float)cfg.P8[PIDALT] * (AltHold - EstAlt) * 0.005f);
    BaroI  = (int16_t)((float)cfg.I8[PIDALT] * vario * 0.02f);   // That is actually a "D"
    BaroD  = (int16_t)((float)cfg.D8[PIDALT] * (100.0f - ThrAngle) * 0.04f); // That is actually the Tiltcompensation
}
#endif

/*
/////// GPS INS TESTCODE
//  Testcode
    static uint32_t previous5HzT = 0;
		flthead = 0;                                                // if no mag do bodyframe below
//  Testcode
    int16_t knob = constrain(rcData[AUX3]-1000,0,1000);
		float  knobbi= (float)knob * 0.001f;
		debug[0] = knobbi * 1000;
	  if (currentT > previous5HzT + 200000){
        previous5HzT = currentT;
		    VelNorth     = VelNorth * knobbi;
        VelEast      = VelEast  * knobbi;
		}
		debug[1] = VelNorth;
		debug[2] = VelEast;
/////// GPS INS TESTCODE
*/
