@echo off
REM ============================================================
REM Squidiff 训练脚本 - Mast Cell 数据 (Windows)
REM ============================================================

REM 设置项目根目录
set PROJECT_ROOT=E:\Development\Squidiff
cd /d %PROJECT_ROOT%

REM ============================================================
REM 参数配置
REM ============================================================

REM 数据路径
set DATA_PATH=data\mast_cells_200genes.h5ad
set GENE_SIZE=200
set OUTPUT_DIM=200

REM 日志路径 (使用时间戳)
for /f "tokens=2-4 delims=/ " %%a in ('date /t') do (set mydate=%%c%%a%%b)
for /f "tokens=1-2 delims=: " %%a in ('time /t') do (set mytime=%%a%%b)
set LOGGER_PATH=logs\mast_cells_%mydate%_%mytime%

REM 模型保存路径
set RESUME_CHECKPOINT=mast_cells_model

REM 训练参数
set BATCH_SIZE=64
set MICROBATCH=-1
set LR=1e-4
set WEIGHT_DECAY=0.0
set LR_ANNEAL_STEPS=50000
set EMA_RATE=0.9999

REM Diffusion 参数
set DIFFUSION_STEPS=1000
set NOISE_SCHEDULE=linear
set LEARN_SIGMA=FALSE

REM 模型参数
set NUM_LAYERS=3
set NUM_CHANNELS=128
set DROPOUT=0.0
set USE_FP16=FALSE
set USE_ENCODER=TRUE

REM 其他参数
set USE_DRUG_STRUCTURE=FALSE
set LOG_INTERVAL=1000
set SAVE_INTERVAL=10000

REM ============================================================
REM 创建日志目录
REM ============================================================
if not exist "%LOGGER_PATH%" mkdir "%LOGGER_PATH%"

REM ============================================================
REM 打印配置
REM ============================================================
echo ============================================================
echo Squidiff 训练配置
echo ============================================================
echo 数据路径:       %DATA_PATH%
echo 基因数:         %GENE_SIZE%
echo 输出维度:       %OUTPUT_DIM%
echo 批大小:         %BATCH_SIZE%
echo 学习率:         %LR%
echo 扩散步数:       %DIFFUSION_STEPS%
echo 日志路径:       %LOGGER_PATH%
echo 检查点:         %RESUME_CHECKPOINT%
echo ============================================================

REM ============================================================
REM 数据验证
REM ============================================================
echo.
echo 检查数据文件...
if not exist "%DATA_PATH%" (
    echo 错误: 数据文件不存在: %DATA_PATH%
    echo 请先运行 scripts\1.data_output.r 生成数据
    pause
    exit /b 1
)
echo √ 数据文件存在

REM ============================================================
REM 开始训练
REM ============================================================
echo.
echo 开始训练...
echo.

python train_squidiff.py ^
  --data_path %DATA_PATH% ^
  --gene_size %GENE_SIZE% ^
  --output_dim %OUTPUT_DIM% ^
  --logger_path %LOGGER_PATH% ^
  --resume_checkpoint %RESUME_CHECKPOINT% ^
  --batch_size %BATCH_SIZE% ^
  --microbatch %MICROBATCH% ^
  --lr %LR% ^
  --weight_decay %WEIGHT_DECAY% ^
  --lr_anneal_steps %LR_ANNEAL_STEPS% ^
  --ema_rate %EMA_RATE% ^
  --diffusion_steps %DIFFUSION_STEPS% ^
  --noise_schedule %NOISE_SCHEDULE% ^
  --learn_sigma %LEARN_SIGMA% ^
  --num_layers %NUM_LAYERS% ^
  --num_channels %NUM_CHANNELS% ^
  --dropout %DROPOUT% ^
  --use_fp16 %USE_FP16% ^
  --use_encoder %USE_ENCODER% ^
  --use_drug_structure %USE_DRUG_STRUCTURE% ^
  --log_interval %LOG_INTERVAL% ^
  --save_interval %SAVE_INTERVAL%

REM ============================================================
REM 训练完成
REM ============================================================
echo.
echo ============================================================
echo √ 训练完成！
echo 模型保存在: %LOGGER_PATH%\%RESUME_CHECKPOINT%
echo ============================================================
pause
