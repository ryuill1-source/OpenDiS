#
# rn matrix
# col(1): pos x, col(2): pos y, col(3): pos z, col(4): flag, 
# col(5): node #, col(6): number of neighbors
#
# links matrix
# col(1): node #, col(2): neighbor node #,
# col(3): Burgers x, col(4): Burgers y, col(5): Burgers z
# col(6): normal x,  col(7): normal y,  col(8): normal z
#
# restart.data
#	Primary lines: 0,node_tag, x, y, z, num_arms, constraint
#	Secondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz
import re

def generate_rnlinks2(filename):
    #filename='restart.data'

    fileID=open(filename,'r')
    keyword="Secondary lines:"
    # 1. Find the line number (num_data_start) data begins -1
    for i,lines in enumerate(fileID):
        if keyword in lines:
            break
    # 2. Generate 'tot_data' dictionart from data file
    tot_data = {}  # Dictionary data with [node_id] : [node_data, seg1_data, seg2_data ... ]
    node_list = []
    for i,lines in enumerate(fileID):
            temp = lines.strip()
            str_list = re.split(r'[,\s]+',temp)   # Using re.split for the case there is whitespace between CPU ID, node Id
            data_list = [float(x) for x in str_list]
            if len(data_list)==0:
                break
            elif len(data_list)==7:
                node_id = data_list[1]
                node_list.append(node_id)
                node_data = data_list[2:]
                tot_data[node_id]=[node_data]
            elif len(data_list)==5:
                seg_data = [node_id]+data_list[1:]
            elif len(data_list)==3:
                seg_data.extend(data_list)
                tot_data[node_id].append(seg_data)
    
    # 3. Generate 'rn matrix' from data file
    rn=[]
    num_node=0
    for node_id in node_list:
        if len(tot_data[node_id])==1:   # broken node
            print(f'    + Broken link detected at node {node_id:.0f} ')
        else:
            node_data=tot_data[node_id][0]
            rn_data=node_data[:3]+[node_data[-1]]+[node_id]+[node_data[-2]]
            rn.append(rn_data)
    real_rn_id=[rn_data[4] for rn_data in rn]

    # 4. Generate 'links matrix' from data file
    links=[]
    for node_id in node_list:
        if len(tot_data[node_id])>1:  
            seg_list=tot_data[node_id][1:]
            for seg in seg_list:
                # 5. Remove repetations in links data.
                if seg[0] < seg[1]:
                    links.append(seg)
                    # 6. Correct links if ghost rn exists
                    if seg[0] not in real_rn_id and seg[1] not in real_rn_id:
                        print(f'    + Link ({seg[0]:.0f},{seg[1]:.0f}) has unknown node.')

    id_1 = [row[0] for row in links]
    id_2 = [row[1] for row in links]
    link_id = id_1+id_2 

    # 7. Check neighbor numbers in rn and correct if it is incorrect.
    for i,rn_data in enumerate(rn):
        node_id=rn_data[-2]
        seg_num=float(link_id.count(node_id))
        if seg_num != rn_data[-1]:
            print(f'    + Arm number {rn_data[-1]:.0f} in node {node_id:.0f} is wrong. ({seg_num:.0f} is correct!)')
            rn[i][-1]=seg_num         
    
    fileID.close()
    return rn, links

