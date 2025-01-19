#!/bin/bash

script_dir=$(dirname "$0")
src_dir="$script_dir/../src"
venv_dir="$script_dir/.venv"

source "$venv_dir/bin/activate"

# Find project directories containing main.cpp
find "$src_dir" -name "main.cpp" -print0 | while IFS= read -r -d $'\0' main_cpp; do
  project_dir=$(dirname "$main_cpp")
  project_name=$(basename "$project_dir")

  echo "Benchmarking project: $project_name"
  python3 "$script_dir/mac_bench.py" "$project_name"
done

deactivate
