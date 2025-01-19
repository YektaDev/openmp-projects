import subprocess
import matplotlib.pyplot as plt
import os

def run_experiment(project_name, build_script, executable_name, num_threads_list, output_file):
    """Runs the experiment for a given build and executable."""
    subprocess.run([build_script, project_name], check=True, shell=True)

    executable_path = os.path.join("../../build", project_name, executable_name)
    for num_threads in num_threads_list:
        subprocess.run([executable_path, str(num_threads), output_file], check=True)

def plot_results(debug_file, release_file):
    """Plots the results from the debug and release data files."""
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

    plt.plot(num_threads_debug, execution_times_debug, marker='o', label='Debug')
    plt.plot(num_threads_release, execution_times_release, marker='x', label='Release')
    plt.xlabel("Number of Threads")
    plt.ylabel("Execution Time (seconds)")
    plt.title("Execution Time vs. Number of Threads (Debug vs. Release)")
    plt.legend()
    plt.grid(True)
    plt.savefig("output.png")
    plt.show()

if __name__ == "__main__":
    project_name = "your_project_name"  # Replace with your project name
    num_threads_list = [1, 2, 4, 8, 16]
    debug_output_file = "debug_results.txt"
    release_output_file = "release_results.txt"

    # Clean previous results (optional)
    if os.path.exists(debug_output_file):
        os.remove(debug_output_file)
    if os.path.exists(release_output_file):
        os.remove(release_output_file)

    run_experiment(project_name, "./mac_debug.sh", "main-debug", num_threads_list, debug_output_file)
    run_experiment(project_name, "./mac_release.sh", "main-release", num_threads_list, release_output_file)

    plot_results(debug_output_file, release_output_file)