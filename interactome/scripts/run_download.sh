# Start with run_sampling
timestamp=$(date "+%Y-%m-%d_%H-%M-%S")
run_dir="$LABHOME/interactome/logs/download/runs/${timestamp}"
mkdir -p "$run_dir"

CUDA_VISIBLE_DEVICES=0 nohup python -u -m scripts.download \
  hydra.run.dir=${run_dir} \
  download="pepmlm" \
  > ${run_dir}/run.log 2>&1 &

echo $! > ${run_dir}/pid.txt
echo "Run started with PID $(cat ${run_dir}/pid.txt)   |   Logs at ${run_dir}/run.log"