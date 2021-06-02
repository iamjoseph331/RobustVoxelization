import os
import sys
import subprocess

filenames = ["100026", "100643", "102041", "105803", "113001", "1207662",  "123044", "125071", "126286", "128915", "1344039", "1344043", "1344053", "1344057", "137977", "138203", "139934", "1439537"]
dir_ftetwilds = []

dir_ftetwilds.append("error/error_1_2e_2")
dir_ftetwilds.append("error/error_1e_2")
dir_ftetwilds.append("error/error_5e_3")
dir_ftetwilds.append("error/error_2e_3")
dir_ftetwilds.append("error/error_1e_3")


for dir_ftetwild in dir_ftetwilds:
    print(dir_ftetwild)
    for filename in filenames:
        command =  ["./exe.out", "-r", "128", "-i", os.path.join(dir_ftetwild, filename+".msh"), "-o", os.path.join(dir_ftetwild, filename), "-p", "-v"]
        #print("Start" + filename)
        subprocess.run(command)
        #print("End" + filename)

dir_tetgen = "error/tetgen_ele_node"
print(dir_tetgen)
for filename in filenames:
    command =  ["./exe.out", "-n", "-r", "128", "-i", os.path.join(dir_tetgen, filename), "-o", os.path.join(dir_tetgen, filename), "-p", "-v"]
    subprocess.run(command)
    print(filename)