import os
import sys

filenames = ["100026", "100643", "102041", "105803", "113001", "1207662",  "123044", "125071", "126286", "128915", "1344039", "1344043", "1344053", "1344057", "137977", "138203", "139934", "1439537"]

dir_tetgen = "error/tetgen_ele_node"
dir_ftetwilds = []

dir_ftetwilds.append("error/error_1_2e_2")
dir_ftetwilds.append("error/error_1e_2")
dir_ftetwilds.append("error/error_5e_3")
dir_ftetwilds.append("error/error_2e_3")
dir_ftetwilds.append("error/error_1e_3")
# errors = ["2e-3", "1e-3", "5e-4", "1e-4"]
errors = ["1.2e-2", "1e-2",  "5e-3", "2e-3", "1e-3"]

result = open('result_128_ver2.csv', 'w')
# result.write('Num,2e-3,1e-3,5e-4,1e-4,\n')
result.write('Num,1.2e-2,1e-2,5e-3,2e-3,1e-3,\n')

# filename = "100026"
# dir_ftetwild = "/home/inm/master/08/r08944020/exp/RobustVoxelization/ftetwild_msh"
for filename in filenames:
    if not os.path.isfile(os.path.join(dir_tetgen, filename)):
        print("Can't find {} in tetgen_ele_node".format(filename))
        continue
    
    print("file number:", filename)
    j = 0
    cnt_diff = [0 for _ in range(5)]
    
    file_tetgen = open(os.path.join(dir_tetgen, filename), 'r')
    file_tetgen.readline()
    nums_tetgen = [x for x in file_tetgen.readline().split(" ")] 
    cnt = 0
    for i in range(len(nums_tetgen)):
        if nums_tetgen[i] == '1':
            cnt += 1
    for dir_ftetwild in dir_ftetwilds:        
        file_ftetwild = open(os.path.join(dir_ftetwild, filename), 'r')    
        file_ftetwild.readline()
        nums_ftetwild = [x for x in file_ftetwild.readline().split(" ")]             

        for i in range(len(nums_tetgen)):
            if nums_tetgen[i] != nums_ftetwild[i]:
                cnt_diff[j] += 1

        print("{}: num of differences = {}".format(errors[j], cnt_diff[j]))
        j += 1
    result.write('{},{},{},{},{},{},{}\n'.format(filename, cnt_diff[0], cnt_diff[1], cnt_diff[2], cnt_diff[3], cnt_diff[4], cnt))
    print('ToTal:' + str(cnt))


