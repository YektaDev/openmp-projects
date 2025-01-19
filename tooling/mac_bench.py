import subprocess
import matplotlib.pyplot as plt
import os
import sys

def run_experiment(project_name, build_script, executable_name, num_threads_list, output_file, script_dir, benchmark_dir):
    """Runs the experiment for a given build and executable."""
    subprocess.run([os.path.join(script_dir, build_script), project_name], check=True, shell=False)

    executable_path = os.path.abspath(os.path.join(script_dir, "..", "build", project_name, executable_name))
    for num_threads in num_threads_list:
        # Redirect output to a file in the benchmark directory
        output_file_path = os.path.join(benchmark_dir, output_file)
        subprocess.run([executable_path, str(num_threads), output_file_path], check=True)

def plot_results(project_name, debug_file, release_file, benchmark_dir):
    """Plots the results and saves the chart as a high-quality PNG."""
    def read_data(data_file):
        if not os.path.exists(data_file):
            print(f"Warning: Data file not found: {data_file}")
            return [], [], [], [], []

        num_threads = []
        avg_times = []
        min_times = []
        max_times = []
        std_devs = []
        with open(data_file, "r") as f:
            for line in f:
                parts = line.split()
                if len(parts) == 5:
                    num_threads.append(int(parts[0]))
                    avg_times.append(float(parts[1]))
                    min_times.append(float(parts[2]))
                    max_times.append(float(parts[3]))
                    std_devs.append(float(parts[4]))
        return num_threads, avg_times, min_times, max_times, std_devs

    num_threads_debug, avg_times_debug, min_times_debug, max_times_debug, std_devs_debug = read_data(debug_file)
    num_threads_release, avg_times_release, min_times_release, max_times_release, std_devs_release = read_data(release_file)

    if not num_threads_debug or not num_threads_release:
        print("Skipping plot due to missing data.")
        return

    # --- Styling for a more beautiful chart ---
    plt.figure(figsize=(12, 7), dpi=300)  # Customize figure size and DPI

    # Plot with error bars (using standard deviation)
    plt.errorbar(num_threads_debug, avg_times_debug, yerr=std_devs_debug, fmt='-o', color='#007acc', label='Debug', capsize=5)
    plt.errorbar(num_threads_release, avg_times_release, yerr=std_devs_release, fmt='-x', color='#ff7f0e', label='Release', capsize=5)

    plt.xlabel("Number of Threads", fontsize=12, fontweight='bold')
    plt.ylabel("Execution Time (seconds)", fontsize=12, fontweight='bold')
    plt.title(f"Execution Time vs. Number of Threads ({project_name})", fontsize=14, fontweight='bold')

    plt.legend(fontsize=10, loc='upper left', frameon=True, shadow=True)
    plt.grid(True, linestyle='--', alpha=0.7)  # Add a subtle grid

    # Customize tick parameters
    plt.tick_params(axis='both', which='major', labelsize=10, direction='inout', length=6, width=1)

    # Add a light background color
    ax = plt.gca()
    ax.set_facecolor('#f8f8f8')

    # Tight layout to minimize whitespace
    plt.tight_layout()

    # --- Save the plot ---
    output_filename = os.path.join(benchmark_dir, f"{project_name}.png")
    plt.savefig(output_filename, bbox_inches='tight')
    plt.close()

    print(f"Benchmark chart saved to: {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python mac_bench.py <project_name>")
        sys.exit(1)

    project_name = sys.argv[1]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    num_threads_list = [1, 2, 4, 8, 12, 16, 20, 24, 28, 32]
    debug_output_file = f"{project_name}_debug.txt"
    release_output_file = f"{project_name}_release.txt"

    # Define the benchmark directory
    benchmark_dir = os.path.abspath(os.path.join(script_dir, "..", "benchmark"))

    # Create the benchmark directory if it doesn't exist
    os.makedirs(benchmark_dir, exist_ok=True)

    # Clean previous results in the benchmark directory
    debug_file_path = os.path.join(benchmark_dir, debug_output_file)
    release_file_path = os.path.join(benchmark_dir, release_output_file)
    if os.path.exists(debug_file_path):
        os.remove(debug_file_path)
    if os.path.exists(release_file_path):
        os.remove(release_file_path)

    run_experiment(project_name, "mac_debug.sh", "main-debug", num_threads_list, debug_output_file, script_dir, benchmark_dir)
    run_experiment(project_name, "mac_release.sh", "main-release", num_threads_list, release_output_file, script_dir, benchmark_dir)

    plot_results(project_name, debug_file_path, release_file_path, benchmark_dir)