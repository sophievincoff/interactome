# Start with run_sampling
timestamp=$(date "+%Y-%m-%d_%H-%M-%S")
name="alphafolddb_query"
run_dir="$LABHOME/interactome/logs/analyze/runs/${name}/${timestamp}"
mkdir -p "$run_dir"

# process.mode="all" \
CUDA_VISIBLE_DEVICES=0,1 nohup python -u -m scripts.analyze \
  hydra.run.dir=${run_dir} \
  analyze=${name} \
  analyze.workers=32 \
  analyze.debug=false \
  analyze.output_path=interactome/data_files/processed/intact/afdb_query/all_uniprots_afdb_stats.tsv\
  > ${run_dir}/run.log 2>&1 &

echo $! > ${run_dir}/pid.txt
echo "Run started with PID $(cat ${run_dir}/pid.txt)   |   Logs at ${run_dir}/run.log"