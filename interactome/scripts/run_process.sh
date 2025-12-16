# Start with run_sampling
timestamp=$(date "+%Y-%m-%d_%H-%M-%S")
name="pepmlm"
run_dir="$LABHOME/interactome/logs/process/runs/${name}/${timestamp}"
mkdir -p "$run_dir"

# process.mode="all" \
CUDA_VISIBLE_DEVICES=0 nohup python -u -m scripts.process \
  hydra.run.dir=${run_dir} \
  process=${name} \
  process.debug=false \
  > ${run_dir}/run.log 2>&1 &

echo $! > ${run_dir}/pid.txt
echo "Run started with PID $(cat ${run_dir}/pid.txt)   |   Logs at ${run_dir}/run.log"