# Start with run_sampling
timestamp=$(date "+%Y-%m-%d_%H-%M-%S")
run_dir="$HOME/interactome/logs/process/runs/${timestamp}"
mkdir -p "$run_dir"

CUDA_VISIBLE_DEVICES=0 nohup python -u -m scripts.process \
  hydra.run.dir=${run_dir} \
  process="biogrid" \
  > ${run_dir}/run.log 2>&1 &

echo $! > ${run_dir}/pid.txt