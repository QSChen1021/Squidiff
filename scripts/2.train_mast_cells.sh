#!/bin/bash
# ============================================================
# Squidiff 训练脚本 - Mast Cell 数据
# ============================================================
#
# 用途: 使用 Mast Cell 数据训练 Squidiff 模型
# 数据: data/mast_cells_200genes.h5ad (200 基因)
#
# ============================================================

# 设置项目根目录
PROJECT_ROOT="E:/Development/Squidiff"
cd "$PROJECT_ROOT"

# ============================================================
# 参数配置
# ============================================================

# 数据路径
DATA_PATH="data/mast_cells_200genes.h5ad"
GENE_SIZE=200
OUTPUT_DIM=200

# 日志路径
LOGGER_PATH="logs/mast_cells_$(date +%Y%m%d_%H%M%S)"

# 模型保存路径
RESUME_CHECKPOINT="mast_cells_model"

# 训练参数
BATCH_SIZE=64
MICROBATCH=-1
LR=1e-4
WEIGHT_DECAY=0.0
LR_ANNEAL_STEPS=50000
EMA_RATE="0.9999"

# Diffusion 参数
DIFFUSION_STEPS=1000
NOISE_SCHEDULE="linear"
LEARN_SIGMA=FALSE

# 模型参数
NUM_LAYERS=3
NUM_CHANNELS=128
DROPOUT=0.0
USE_FP16=FALSE
USE_ENCODER=TRUE

# 其他参数
USE_DRUG_STRUCTURE=FALSE
LOG_INTERVAL=1000
SAVE_INTERVAL=10000

# ============================================================
# 创建日志目录
# ============================================================
mkdir -p "$LOGGER_PATH"

# ============================================================
# 打印配置
# ============================================================
echo "============================================================"
echo "Squidiff 训练配置"
echo "============================================================"
echo "数据路径:       $DATA_PATH"
echo "基因数:         $GENE_SIZE"
echo "输出维度:       $OUTPUT_DIM"
echo "批大小:         $BATCH_SIZE"
echo "学习率:         $LR"
echo "扩散步数:       $DIFFUSION_STEPS"
echo "日志路径:       $LOGGER_PATH"
echo "检查点:         $RESUME_CHECKPOINT"
echo "============================================================"

# ============================================================
# 数据验证
# ============================================================
echo ""
echo "检查数据文件..."
if [ ! -f "$DATA_PATH" ]; then
    echo "错误: 数据文件不存在: $DATA_PATH"
    echo "请先运行 scripts/1.data_output.r 生成数据"
    exit 1
fi
echo "✓ 数据文件存在"

# ============================================================
# 开始训练
# ============================================================
echo ""
echo "开始训练..."
echo ""

python train_squidiff.py \
  --data_path "$DATA_PATH" \
  --gene_size $GENE_SIZE \
  --output_dim $OUTPUT_DIM \
  --logger_path "$LOGGER_PATH" \
  --resume_checkpoint "$RESUME_CHECKPOINT" \
  --batch_size $BATCH_SIZE \
  --microbatch $MICROBATCH \
  --lr $LR \
  --weight_decay $WEIGHT_DECAY \
  --lr_anneal_steps $LR_ANNEAL_STEPS \
  --ema_rate $EMA_RATE \
  --diffusion_steps $DIFFUSION_STEPS \
  --noise_schedule $NOISE_SCHEDULE \
  --learn_sigma $LEARN_SIGMA \
  --num_layers $NUM_LAYERS \
  --num_channels $NUM_CHANNELS \
  --dropout $DROPOUT \
  --use_fp16 $USE_FP16 \
  --use_encoder $USE_ENCODER \
  --use_drug_structure $USE_DRUG_STRUCTURE \
  --log_interval $LOG_INTERVAL \
  --save_interval $SAVE_INTERVAL

# ============================================================
# 训练完成
# ============================================================
echo ""
echo "============================================================"
echo "✓ 训练完成！"
echo "模型保存在: $LOGGER_PATH/$RESUME_CHECKPOINT"
echo "============================================================"
