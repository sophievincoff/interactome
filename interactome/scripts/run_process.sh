# Start with run_sampling
timestamp=$(date "+%Y-%m-%d_%H-%M-%S")
name="pinderdb"
run_dir="$HOME/interactome/logs/process/runs/${name}/${timestamp}"
mkdir -p "$run_dir"

CUDA_VISIBLE_DEVICES=0,1,2,3 nohup python -u -m scripts.process \
  hydra.run.dir=${run_dir} \
  process=${name} \
  process.mode="positive" \
  process.n_cores=22 \
  process.per_core_samples=30000 \
  process.max_workers_inner=1 \
  process.cover_all=true \
  > ${run_dir}/run.log 2>&1 &

echo $! > ${run_dir}/pid.txt