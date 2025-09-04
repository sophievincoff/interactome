# Start with run_sampling
timestamp=$(date "+%Y-%m-%d_%H-%M-%S")
run_dir="$HOME/interactome/logs/download/runs/${timestamp}"
mkdir -p "$run_dir"

CUDA_VISIBLE_DEVICES=0 nohup python -u -m scripts.download \
  hydra.run.dir=${run_dir} \
  download="biogrid" \
  > ${run_dir}/run.log 2>&1 &

echo $! > ${run_dir}/pid.txt