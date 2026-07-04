import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt

# Get current directory
current_dir = os.getcwd()
Vg_max = 3
J = []
Vg = []
n =50
for i in range(n+1):
    datavec = []
    V_g = i*Vg_max/50
    Vg.append(V_g)
    with open("./data","r") as f:
        for line in f:
            datavec.append(line)
    datavec[3] = str(V_g)+"\n"
    with open("./data","w") as f:
        f.writelines(datavec)
    with open("./Vg.txt","w") as g:
        g.write(str(V_g))
# List of C files to compile and run
    c_files = ["poisson.c","T_E.c","gatecurrent.c"]  # replace with your actual filenames

    for c_file in c_files:
        if not os.path.exists(os.path.join(current_dir, c_file)):
            print(f"❌ {c_file} not found in {current_dir}")
            continue

        exe_file = c_file.replace(".c", "")

        print(f"\n🔧 Compiling {c_file}...")
        compile_cmd = ["gcc", c_file, "-o", exe_file,"-lm"]
        compile_result = subprocess.run(compile_cmd, capture_output=True, text=True)

        if compile_result.returncode != 0:
            print(f"❌ Compilation failed for {c_file}:")
            print(compile_result.stderr)
            continue
        else:
            print(f"✅ Compilation successful: {exe_file}")

        print(f"🚀 Running {exe_file}...")
        run_result = subprocess.run([f"./{exe_file}"], capture_output=True, text=True)
    # if (c_files == c_files[2]):
        if run_result.returncode == 0:
            print(f"✅ Output of {exe_file}:\n{run_result.stdout}")
        else:
            print(f"⚠️ Runtime error in {exe_file}:\n{run_result.stderr}")
        file = open("./output.txt", "r")
        content = -1*np.double(file.read().strip())
        print(content)
        if(c_file == c_files[2]):
            J.append(content)
        file.close()

print(Vg)
print(J)
plt.scatter(Vg,J)
plt.xlabel("Vg(in V)")
plt.ylabel("J(in Acm-2)")
print("\n\n\nrefer to J_Vs_Vg.png for the plot\n")
plt.savefig("J_Vs_Vg.png")
plt.show()
