import subprocess
import matplotlib.pyplot as plt
import os
import sys

def run_experiment(project_name, build_script, executable_name, num_threads_list, output_file, script_dir):
    """Runs the experiment for a given build and executable."""
    subprocess.run([os.path.join(script_dir, build_script), project_name], check=True, shell=False)

    executable_path = os.path.abspath(os.path.join(script_dir, "..", "build", project_name, executable_name))
    for num_threads in num_threads_list:
        subprocess.run([executable_path, str(num_threads), output_file], check=True)

def plot_results(project_name, debug_file, release_file, benchmark_dir):
    """Plots the results and saves the chart as a high-quality PNG."""
    def read_data(data_file):
        num_threads = []
        execution_times = []
        with open(data_file, "r") as f:
            for line in f:
                parts = line.split()
                if len(parts) == 2:
                    num_threads.append(int(parts[0]))
                    execution_times.append(float(parts[1]))
        return num_threads, execution_times

    num_threads_debug, execution_times_debug = read_data(debug_file)
    num_threads_release, execution_times_release = read_data(release_file)

    plt.figure(figsize=(10, 6), dpi=300)  # Customize figure size and DPI
    plt.plot(num_threads_debug, execution_times_debug, marker='o', label='Debug')
    plt.plot(num_threads_release, execution_times_release, marker='x', label='Release')
    plt.xlabel("Number of Threads")
    plt.ylabel("Execution Time (seconds)")
    plt.title(f"Execution Time vs. Number of Threads (Debug vs. Release) - {project_name}")
    plt.legend()
    plt.grid(True)

    # Save the plot to the benchmark directory
    output_filename = os.path.join(benchmark_dir, f"{project_name}_benchmark.png")
    plt.savefig(output_filename, bbox_inches='tight')  # Adjust layout
    plt.close()  # Close the figure to avoid displaying it

    print(f"Benchmark chart saved to: {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python mac_bench.py <project_name>")
        sys.exit(1)

    project_name = sys.argv[1]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    num_threads_list = [1, 2, 4, 8, 16]
    debug_output_file = f"{project_name}_debug_results.txt"
    release_output_file = f"{project_name}_release_results.txt"

    # Define the benchmark directory
    benchmark_dir = os.path.abspath(os.path.join(script_dir, "..", "benchmark"))

    # Create the benchmark directory if it doesn't exist
    os.makedirs(benchmark_dir, exist_ok=True)

    # Clean previous results (optional)
    if os.path.exists(debug_output_file):
        os.remove(debug_output_file)
    if os.path.exists(release_output_file):
        os.remove(release_output_file)

    run_experiment(project_name, "mac_debug.sh", "main-debug", num_threads_list, debug_output_file, script_dir)
    run_experiment(project_name, "mac_release.sh", "main-release", num_threads_list, release_output_file, script_dir)

    plot_results(project_name, debug_output_file, release_output_file, benchmark_dir)
