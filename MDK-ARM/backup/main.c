/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2026 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "dac.h"
#include "i2c.h"
#include "tim.h"
#include "gpio.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include <stdio.h>
#include <string.h>
#include "lsm6dsv_reg.h"
#include "stdbool.h"
#include "math.h"
#include "vqf6d.h"
#include "Estimate.h"
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */
#define DAC_MAX_VAL 0xFFF // 
#define I2C1_SCL_PORT          GPIOB
#define I2C1_SCL_PIN           GPIO_PIN_6

#define I2C1_SDA_PORT          GPIOB
#define I2C1_SDA_PIN           GPIO_PIN_7			

#define I2C2_SCL_PORT          GPIOB
#define I2C2_SCL_PIN           GPIO_PIN_10

#define I2C2_SDA_PORT          GPIOB
#define I2C2_SDA_PIN           GPIO_PIN_11	

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */
int32_t writeI2CReg(void *handle, uint8_t reg, const uint8_t *bufp, uint16_t len);
int32_t readI2CReg(void *handle, uint8_t reg, uint8_t *bufp, uint16_t len);
void delay_ms(uint32_t ms);
void GetRawSensor(stmdev_ctx_t *dev_ctx,int16_t *acc_raw, int16_t *gyro_raw);
void GetCorrectedSensor(stmdev_ctx_t *dev_ctx, float *corrected_acc, float *corrected_gyr, int16_t *gyr_offset);
void ConfigureIMU(void);
HAL_StatusTypeDef I2C_BusRecovery(I2C_HandleTypeDef *I2Cx);
void CalibGyro(void);
/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */
uint8_t whoAmI1 = 0xFF;
uint8_t whoAmI2 = 0xFF;
uint8_t rst = 0;
uint8_t cnt1s = 0;
static int32_t acc_raw_acum1[3],acc_raw_acum2[3];
static int32_t gyr_raw_acum1[3],gyr_raw_acum2[3];
int16_t gyr_offset1[3] = {0,0,0};
int16_t gyr_offset2[3] = {0,0,0};
volatile int16_t acc_offset[3] = {0, 0, 0};
uint16_t calib_step = 0;
float compute_time;
//int32_t test_sum_gyr[3] = {0, 0, 0};
//uint16_t test_count = 0;
//float test_bias_x,test_bias_y,test_bias_z;
bool calib_enable = false; bool calib_occur = false;
float gValue1[3],gValue2[3];
float gyro1[3],gyro2[3];
float acc_test = 0;
float dac_bending = 0;
float dac_twisting = 0;
float test_compute1 = 0;
float test_compute2 = 0;
HAL_StatusTypeDef status;
uint16_t i2c2_error = 0;
float q1_E1S1[4], q2_E2S2[4];
float R1_E1S1[3][3], R2_E2S2[3][3];
static lsm6dsv_filt_settling_mask_t filt_settling_mask;
lsm6dsv_all_sources_t all_sources1,all_sources2;
stmdev_ctx_t dev_ctx1, dev_ctx2;
HAL_StatusTypeDef status_i2c2;
uint32_t error_i2c2;
static VQF6D vqf1;
static VQF6D vqf2;
Estimator estimator;

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */
/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
	

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DAC_Init();
  MX_I2C1_Init();
  MX_TIM1_Init();
  MX_TIM2_Init();
  MX_I2C2_Init();
  MX_TIM3_Init();
  MX_TIM4_Init();
  /* USER CODE BEGIN 2 */
	HAL_TIM_Base_Start(&htim1);
	HAL_TIM_Base_Start_IT(&htim2);
	HAL_TIM_Base_Start(&htim4);
	HAL_DAC_Start(&hdac, DAC_CHANNEL_1);
	dev_ctx1.write_reg = writeI2CReg;
	dev_ctx1.read_reg = readI2CReg;
	dev_ctx1.handle = &hi2c1;
	dev_ctx1.mdelay = delay_ms;
	dev_ctx2.write_reg = writeI2CReg;
	dev_ctx2.read_reg = readI2CReg;
	dev_ctx2.handle = &hi2c2;
	dev_ctx2.mdelay = delay_ms;
//	uint8_t i;
//	for(i = 1; i < 128; i++)
//	{
//		 if(HAL_I2C_IsDeviceReady(&hi2c2, (i << 1), 1, 10) == HAL_OK) 
//		 {
//			 test_compute1 = 99;
//			 test_compute2 = i;
//			 break;
//     }
//		 else
//		 {
//				i2c2_error = HAL_I2C_GetError(&hi2c2);
//		 }
//	}
	lsm6dsv_device_id_get(&dev_ctx1, &whoAmI1);
 	while(whoAmI1 != LSM6DSV_ID)
	{
		rst = I2C_BusRecovery(&hi2c1);
		lsm6dsv_device_id_get(&dev_ctx1, &whoAmI1);
	}
	rst = lsm6dsv_sw_reset(&dev_ctx1);
	delay_ms(10);
	lsm6dsv_device_id_get(&dev_ctx2, &whoAmI2);
 	while(whoAmI2 != LSM6DSV_ID)
	{
		rst = I2C_BusRecovery(&hi2c2);
		lsm6dsv_device_id_get(&dev_ctx2, &whoAmI2);
	}
	rst = lsm6dsv_sw_reset(&dev_ctx2);
	delay_ms(10);
	/* Restore default configuration */
	lsm6dsv_sw_por(&dev_ctx1);
	delay_ms(10);
	lsm6dsv_sw_por(&dev_ctx2);
	delay_ms(10);
	ConfigureIMU();
	delay_ms(50);
	CalibGyro();
	VQF6D_Params p;
	VQF6D_GetDefaultParams(&p);
	VQF6D_Init(&vqf1, 0.01f, 0.01f, &p);
	VQF6D_Init(&vqf2, 0.01f, 0.01f, &p);
	InitEstimator(&estimator);
	GetCorrectedSensor(&dev_ctx1,gValue1,gyro1, gyr_offset1);
	GetCorrectedSensor(&dev_ctx2,gValue2,gyro2, gyr_offset2);
	// Run VQF Sensor Fushion
	VQF6D_Update(&vqf1, gyro1, gValue1);
	VQF6D_GetQuat6D(&vqf1, q1_E1S1);
	VQF6D_QuatToRotMat(q1_E1S1, R1_E1S1);
	VQF6D_Update(&vqf2, gyro2, gValue2);
	VQF6D_GetQuat6D(&vqf2, q2_E2S2);
	VQF6D_QuatToRotMat(q2_E2S2, R2_E2S2);
	memcpy(estimator.debug.R1_ref, R1_E1S1, sizeof(estimator.debug.R1_ref));
	memcpy(estimator.debug.R2_ref, R2_E2S2, sizeof(estimator.debug.R2_ref));
	HAL_TIM_Base_Start_IT(&htim3);
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE1);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLM = 8;
  RCC_OscInitStruct.PLL.PLLN = 168;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLQ = 4;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV4;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV2;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_5) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */
int32_t writeI2CReg(void *handle, uint8_t reg, const uint8_t *bufp, uint16_t len)
{
	HAL_I2C_Mem_Write(handle, LSM6DSV_I2C_ADD_H, reg, I2C_MEMADD_SIZE_8BIT, (uint8_t*) bufp, len, 1000);
	return 0;
}
int32_t readI2CReg(void *handle, uint8_t reg, uint8_t *bufp, uint16_t len)
{
	HAL_I2C_Mem_Read(handle, LSM6DSV_I2C_ADD_H, reg, I2C_MEMADD_SIZE_8BIT, bufp, len, 1000);
	return 0;
}
void delay_ms(uint32_t ms)
{
	__HAL_TIM_SET_COUNTER(&htim1,0); 
	while (__HAL_TIM_GET_COUNTER(&htim1) < ms*10);
}

void ConfigureIMU(void)
{
    /* Enable Block Data Update */
    lsm6dsv_block_data_update_set(&dev_ctx1, PROPERTY_ENABLE);
    lsm6dsv_xl_data_rate_set(&dev_ctx1, LSM6DSV_ODR_AT_1920Hz);
    lsm6dsv_gy_data_rate_set(&dev_ctx1, LSM6DSV_ODR_AT_1920Hz);

    /* Set full scale */
    lsm6dsv_xl_full_scale_set(&dev_ctx1, LSM6DSV_4g);
    lsm6dsv_gy_full_scale_set(&dev_ctx1, LSM6DSV_500dps);
    /*Set mode*/
    lsm6dsv_xl_mode_set(&dev_ctx1, LSM6DSV_XL_HIGH_PERFORMANCE_MD);
    lsm6dsv_gy_mode_set(&dev_ctx1, LSM6DSV_GY_HIGH_PERFORMANCE_MD);
    /* Configure filtering chain */
    filt_settling_mask.drdy = PROPERTY_ENABLE;
    filt_settling_mask.irq_xl = PROPERTY_ENABLE;
    filt_settling_mask.irq_g = PROPERTY_ENABLE;
    lsm6dsv_filt_settling_mask_set(&dev_ctx1, filt_settling_mask);
	// Gyrocope LPF Bandwidth: 154Hz
    lsm6dsv_filt_gy_lp1_set(&dev_ctx1, PROPERTY_ENABLE);
    lsm6dsv_filt_gy_lp1_bandwidth_set(&dev_ctx1, LSM6DSV_GY_VERY_LIGHT);
	// Accelerometer LPF Bandwidth: 
    lsm6dsv_filt_xl_lp2_set(&dev_ctx1, PROPERTY_ENABLE);
    lsm6dsv_filt_xl_lp2_bandwidth_set(&dev_ctx1, LSM6DSV_XL_VERY_LIGHT);
		
		 /* Enable Block Data Update */
    lsm6dsv_block_data_update_set(&dev_ctx2, PROPERTY_ENABLE);
    lsm6dsv_xl_data_rate_set(&dev_ctx2, LSM6DSV_ODR_AT_1920Hz);
    lsm6dsv_gy_data_rate_set(&dev_ctx2, LSM6DSV_ODR_AT_1920Hz);

    /* Set full scale */
    lsm6dsv_xl_full_scale_set(&dev_ctx2, LSM6DSV_4g);
    lsm6dsv_gy_full_scale_set(&dev_ctx2, LSM6DSV_500dps);
    /*Set mode*/
    lsm6dsv_xl_mode_set(&dev_ctx2, LSM6DSV_XL_HIGH_PERFORMANCE_MD);
    lsm6dsv_gy_mode_set(&dev_ctx2, LSM6DSV_GY_HIGH_PERFORMANCE_MD);
    lsm6dsv_filt_settling_mask_set(&dev_ctx2, filt_settling_mask);
	// Gyrocope LPF Bandwidth: 154Hz
    lsm6dsv_filt_gy_lp1_set(&dev_ctx2, PROPERTY_ENABLE);
    lsm6dsv_filt_gy_lp1_bandwidth_set(&dev_ctx2, LSM6DSV_GY_MEDIUM);
	// Accelerometer LPF Bandwidth: 
    lsm6dsv_filt_xl_lp2_set(&dev_ctx2, PROPERTY_ENABLE);
    lsm6dsv_filt_xl_lp2_bandwidth_set(&dev_ctx2, LSM6DSV_XL_MEDIUM);
}

void GetRawSensor(stmdev_ctx_t *dev_ctx, int16_t *acc_raw, int16_t *gyro_raw)
{
    uint8_t buf[12];
    HAL_StatusTypeDef status;

    status = HAL_I2C_Mem_Read(dev_ctx->handle, LSM6DSV_I2C_ADD_H, LSM6DSV_OUTX_L_G, I2C_MEMADD_SIZE_8BIT, buf, 12, 10);

    if(status != HAL_OK)
    {
				HAL_GPIO_TogglePin(GPIOD, GPIO_PIN_12);
        return;
    }

    gyro_raw[0] = (int16_t)((uint16_t)buf[1] << 8 | buf[0]);
    gyro_raw[1] = (int16_t)((uint16_t)buf[3] << 8 | buf[2]);
    gyro_raw[2] = (int16_t)((uint16_t)buf[5] << 8 | buf[4]);

    acc_raw[0] = (int16_t)((uint16_t)buf[7] << 8 | buf[6]);
    acc_raw[1] = (int16_t)((uint16_t)buf[9] << 8 | buf[8]);
    acc_raw[2] = (int16_t)((uint16_t)buf[11] << 8 | buf[10]);
}

void GetCorrectedSensor(stmdev_ctx_t *dev_ctx, float *corrected_acc, float *corrected_gyr, int16_t *gyr_offset)
{
		int16_t data_raw_acceleration[3] = {0, 0, 0};
		int16_t data_raw_angular_rate[3] = {0, 0, 0};
    GetRawSensor(dev_ctx, data_raw_acceleration, data_raw_angular_rate);
    corrected_acc[0] = lsm6dsv_from_fs4_to_mg(data_raw_acceleration[0])*0.00980665f;
    corrected_acc[1] = lsm6dsv_from_fs4_to_mg(data_raw_acceleration[1])*0.00980665f;
    corrected_acc[2] = lsm6dsv_from_fs4_to_mg(data_raw_acceleration[2])*0.00980665f;
	

    corrected_gyr[0] = lsm6dsv_from_fs500_to_mdps(data_raw_angular_rate[0] - gyr_offset[0])*0.001f ;
    corrected_gyr[1] = lsm6dsv_from_fs500_to_mdps(data_raw_angular_rate[1] - gyr_offset[1])*0.001f ;
    corrected_gyr[2] = lsm6dsv_from_fs500_to_mdps(data_raw_angular_rate[2] - gyr_offset[2])*0.001f ;
	
		corrected_gyr[0] *= 0.0174533;
		corrected_gyr[1] *= 0.0174533;
		corrected_gyr[2] *= 0.0174533;
}
void HAL_TIM_PeriodElapsedCallback(TIM_HandleTypeDef *htim)
{
//    if (htim->Instance == TIM2) 
//    {
////			GetCorrectedSensor(gValue,gyro);
//			dac_bending = (uint16_t)(estimator.bending_angle*22.0611111 + 62);
//			dac_bending = (dac_bending > 4095)? 4095 : (dac_bending < 0) ? 0: dac_bending;
//			dac_twisting = (uint16_t)(estimator.twist_angle*22.0611111 + 62);
//			dac_twisting = (dac_twisting > 4095)? 4095 : (dac_twisting < 0) ? 0: dac_twisting;
//			HAL_DAC_SetValue(&hdac,DAC_CHANNEL_1,DAC_ALIGN_12B_R, dac_bending);
//			HAL_DAC_SetValue(&hdac,DAC_CHANNEL_2,DAC_ALIGN_12B_R, dac_twisting);
//			cnt1s++;
//			if(cnt1s > 100)
//			{
//				HAL_GPIO_TogglePin(GPIOD, GPIO_PIN_12);
//				cnt1s = 0;
//			}
//    }
			if (htim->Instance == TIM3) 
			{
				__HAL_TIM_SET_COUNTER(&htim4, 0);
				//Get Data
				GetCorrectedSensor(&dev_ctx1,gValue1,gyro1, gyr_offset1);
				GetCorrectedSensor(&dev_ctx2,gValue2,gyro2, gyr_offset2);
				// Run VQF Sensor Fushion
				VQF6D_Update(&vqf1, gyro1, gValue1);
				VQF6D_GetQuat6D(&vqf1, q1_E1S1);
				VQF6D_QuatToRotMat(q1_E1S1, R1_E1S1);
				VQF6D_Update(&vqf2, gyro2, gValue2);
				VQF6D_GetQuat6D(&vqf2, q2_E2S2);
				VQF6D_QuatToRotMat(q2_E2S2, R2_E2S2);
				float bias1[3];
				float bias2[3];
				float biasSigma1;
				float biasSigma2;

				float gyro1_corr[3];
				float gyro2_corr[3];

				VQF6D_GetBiasEstimate(&vqf1, bias1, &biasSigma1);
				VQF6D_GetBiasEstimate(&vqf2, bias2, &biasSigma2);

				gyro1_corr[0] = gyro1[0] - bias1[0];
				gyro1_corr[1] = gyro1[1] - bias1[1];
				gyro1_corr[2] = gyro1[2] - bias1[2];

				gyro2_corr[0] = gyro2[0] - bias2[0];
				gyro2_corr[1] = gyro2[1] - bias2[1];
				gyro2_corr[2] = gyro2[2] - bias2[2];
				
				EstimatorDebugUpdate(&estimator, &vqf1, gyro1_corr, gyro2_corr, R1_E1S1, R2_E2S2);

				// Compute J1 Axis every 20Hz
				UpdateJointAxis1(&estimator, gyro1_corr);
				if (!estimator.u2_locked)
				{
					UpdateJointAxis2(&estimator, gyro1_corr, gyro2_corr);
				}
				UpdateHeadingData(&estimator, gyro1_corr, gyro2_corr, R1_E1S1, R2_E2S2);
				
				estimator.compute_heading_counter++;
				if (estimator.compute_heading_counter >= 10U)
				{
					if (!estimator.heading_valid)
					{
							EstimateInitialHeading(&estimator);
					}
					else
					{
							TrackHeadingDelta(&estimator);
					}
					estimator.compute_heading_counter = 0;
					ComputeHeadingDiagnostic(&estimator, gyro1_corr, gyro2_corr, R1_E1S1, R2_E2S2);
				}
//			test_compute1 = atan2f(estimator.ee_omega_meas[0],estimator.ee_omega_meas[1]) * 57.2957795f;
//			test_compute2 = cosf(estimator.q2_pred)*estimator.ee_omega_meas[0] - sinf(estimator.q2_pred)*estimator.ee_omega_meas[1];
				acc_test = sqrt(gValue2[0]*gValue2[0] + gValue2[1]*gValue2[1] + gValue2[2]*gValue2[2]);
				uint16_t dt = __HAL_TIM_GET_COUNTER(&htim4);
				compute_time = dt*0.001f;
			}
}

void CalibGyro(void)
{
		calib_step = 0;
		calib_occur = true;
		gyr_raw_acum1[0] = 0; gyr_raw_acum1[1] = 0; gyr_raw_acum1[2] = 0;
		acc_raw_acum1[0] = 0; acc_raw_acum1[1] = 0; acc_raw_acum1[2] = 0;
	
		gyr_raw_acum2[0] = 0; gyr_raw_acum2[1] = 0; gyr_raw_acum2[2] = 0;
		acc_raw_acum2[0] = 0; acc_raw_acum2[1] = 0; acc_raw_acum2[2] = 0;
		while(calib_occur)
		{
				if(calib_step < 1000)
				{
						__HAL_TIM_SET_COUNTER(&htim4, 0);
						int16_t data_raw_accelerometer[3], data_raw_angular_rate[3];
					  GetRawSensor(&dev_ctx1, data_raw_accelerometer, data_raw_angular_rate);
						gyr_raw_acum1[0] += data_raw_angular_rate[0];
						gyr_raw_acum1[1] += data_raw_angular_rate[1];
						gyr_raw_acum1[2] += data_raw_angular_rate[2];
					
						acc_raw_acum1[0] += data_raw_accelerometer[0];
						acc_raw_acum1[1] += data_raw_accelerometer[1];
						acc_raw_acum1[2] += data_raw_accelerometer[2];
	
						GetRawSensor(&dev_ctx2, data_raw_accelerometer, data_raw_angular_rate);
						gyr_raw_acum2[0] += data_raw_angular_rate[0];
						gyr_raw_acum2[1] += data_raw_angular_rate[1];
						gyr_raw_acum2[2] += data_raw_angular_rate[2];
					
						acc_raw_acum2[0] += data_raw_accelerometer[0];
						acc_raw_acum2[1] += data_raw_accelerometer[1];
						acc_raw_acum2[2] += data_raw_accelerometer[2];
					
						calib_step++;
						
				}
				else
				{
						gyr_offset1[0] = (float)(gyr_raw_acum1[0])/(float)calib_step;
						gyr_offset1[1] = (float)(gyr_raw_acum1[1])/(float)calib_step;
						gyr_offset1[2] = (float)(gyr_raw_acum1[2])/(float)calib_step;
					
						gyr_offset2[0] = (float)(gyr_raw_acum2[0])/(float)calib_step;
						gyr_offset2[1] = (float)(gyr_raw_acum2[1])/(float)calib_step;
						gyr_offset2[2] = (float)(gyr_raw_acum2[2])/(float)calib_step;
					
//						acc_offset[0] = (int16_t)((float)(acc_raw_acum2[0])/calib_step);
//						acc_offset[1] = (int16_t)((float)(acc_raw_acum2[1])/calib_step);
//						acc_offset[2] = (int16_t)((float)(acc_raw_acum2[2])/calib_step);	
						calib_occur = false;
				}
				uint16_t dt = __HAL_TIM_GET_COUNTER(&htim4);
				compute_time = dt*0.001f;
				delay_ms(2);
		}
}
HAL_StatusTypeDef I2C_BusRecovery(I2C_HandleTypeDef *I2Cx)
{
    GPIO_InitTypeDef GPIO_InitStruct = {0};
    
    GPIO_TypeDef *scl_port, *sda_port;
    uint16_t scl_pin, sda_pin;

    if (I2Cx->Instance == I2C1)
    {
        scl_port = I2C1_SCL_PORT;
        scl_pin  = I2C1_SCL_PIN;
        sda_port = I2C1_SDA_PORT;
        sda_pin  = I2C1_SDA_PIN;
    }
    else if (I2Cx->Instance == I2C2)
    {
        scl_port = I2C2_SCL_PORT;
        scl_pin  = I2C2_SCL_PIN;
        sda_port = I2C2_SDA_PORT;
        sda_pin  = I2C2_SDA_PIN;
    }
    else
    {
        return HAL_ERROR;
    }

    HAL_I2C_DeInit(I2Cx);

    GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_OD;
    GPIO_InitStruct.Pull = GPIO_NOPULL;
    GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;

    GPIO_InitStruct.Pin = scl_pin;
    HAL_GPIO_Init(scl_port, &GPIO_InitStruct);

    GPIO_InitStruct.Pin = sda_pin;
    HAL_GPIO_Init(sda_port, &GPIO_InitStruct);

    HAL_GPIO_WritePin(scl_port, scl_pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(sda_port, sda_pin, GPIO_PIN_SET);
    HAL_Delay(1);

    if (HAL_GPIO_ReadPin(scl_port, scl_pin) == GPIO_PIN_RESET)
    {
        goto recovery_failed;
    }

    for (int i = 0; i < 9 && HAL_GPIO_ReadPin(sda_port, sda_pin) == GPIO_PIN_RESET; i++)
    {
        HAL_GPIO_WritePin(scl_port, scl_pin, GPIO_PIN_RESET);
        HAL_Delay(1);

        HAL_GPIO_WritePin(scl_port, scl_pin, GPIO_PIN_SET);
        HAL_Delay(1);

        if (HAL_GPIO_ReadPin(scl_port, scl_pin) == GPIO_PIN_RESET)
        {
            goto recovery_failed;
        }
    }

    HAL_GPIO_WritePin(sda_port, sda_pin, GPIO_PIN_RESET);
    HAL_Delay(1);

    HAL_GPIO_WritePin(scl_port, scl_pin, GPIO_PIN_SET);
    HAL_Delay(1);

    if (HAL_GPIO_ReadPin(scl_port, scl_pin) == GPIO_PIN_RESET)
    {
        goto recovery_failed;
    }

    HAL_GPIO_WritePin(sda_port, sda_pin, GPIO_PIN_SET);
    HAL_Delay(1);

    if (HAL_GPIO_ReadPin(sda_port, sda_pin) != GPIO_PIN_SET ||
        HAL_GPIO_ReadPin(scl_port, scl_pin) != GPIO_PIN_SET)
    {
        goto recovery_failed;
    }

    if (I2Cx->Instance == I2C1) {
        MX_I2C1_Init();
    } else if (I2Cx->Instance == I2C2) {
        MX_I2C2_Init();
    }
    
    return HAL_OK;

recovery_failed:
    if (I2Cx->Instance == I2C1) {
        MX_I2C1_Init();
    } else if (I2Cx->Instance == I2C2) {
        MX_I2C2_Init();
    }
    return HAL_ERROR;
}
/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}
#ifdef USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
