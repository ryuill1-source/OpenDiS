# V14: modified version for indentation
# Not all FEM elements use dislocation plasticity. The stress output from FEM are handled differently
# The Stress_int_point array has 'N_total_elements' rows, but only rowID matching DDM elements has stress info.
import numpy as np
import os
import shutil
from generate_rnlinks2 import generate_rnlinks2
from write_node_data4 import write_node_data4

def generate_postprocessing_data14(Num_Elements,Num_DDM_Elements,Element_type,Elements,Elements_Center,Elements_Min,Elements_Max,
    Jobname,filename_FEMoutput_stress,filename_FEMoutput_nodes,foldername_FEMstress,filename_FEMstress_head,filename_NUCsite_head,foldername_restart,filename_center_element,
    systemtype_Linux,foldername_root,foldername_tests,foldername_temp_win,foldername_temp_wsl,foldername_ABAQUS,EffSlip,dt_FEM_each_step,filename_rn_check,filename_link_element,
    num_cpus,xyzlimit,SurfaceElementSet,Num_SurfaceElement,FEM_Stress_Scale,n_step,enabled_temp_disp):
    
    ## Find nucleation sites based on maximum RSS and maximum distance from the existing dislocation segments
    filename_stress_nucleation='SurfaceStress';     # filename for ParaDiS nucleation
    ## Import stress components for integration points
    Stress_int_point=np.zeros((Num_Elements,8))
    
    for i in range(num_cpus):
        file_stress='./'+foldername_ABAQUS+'/'+filename_FEMoutput_stress+f'{i+1:04d}'+'.txt'
        Stress_int_point_temp=np.loadtxt(file_stress)
        if i==0:
            Stress_int_DDM=Stress_int_point_temp
        else:
            #Stress_int_DDM.extend(Stress_int_point_temp)
            #(25/07/01, seongjin) np.ndarray doesn't have 'extend' 
            Stress_int_DDM = np.vstack([Stress_int_DDM, Stress_int_point_temp])
    
    double_write=0
    if len(Stress_int_DDM)>Num_DDM_Elements:
        double_write=1
    
    for i in range(num_cpus):
        file_stress='./'+foldername_ABAQUS+'/'+filename_FEMoutput_stress+f'{i+1:04d}'+'.txt'
        Stress_int_point_temp=np.loadtxt(file_stress)
        
        if double_write:       # The output from ABAQUS is doubly written in this cycle
            #size_stress=len(Stress_int_point_temp)/2
            #(25/07/01 seongjin) slicing number should be integer type
            size_stress=int(len(Stress_int_point_temp)/2)

            #Stress_int_point_temp[:size_stress]
            #(25/07/01 seongjin) 2D slicing, not 1D
            Stress_int_point_temp = Stress_int_point_temp[:size_stress, :]
        
        if i==0:
            Stress_int_DDM=Stress_int_point_temp
        else:
            #Stress_int_DDM.extend(Stress_int_point_temp)
            #(25/07/01, seongjin) np.ndarray doesn't have 'extend' 
            Stress_int_DDM = np.vstack([Stress_int_DDM, Stress_int_point_temp])

        fileID=open(file_stress,'w')    # Deleting contents
        fileID.close()     
    Stress_int_DDM=np.array(Stress_int_DDM)    
    Stress_int_DDM=Stress_int_DDM[np.argsort(Stress_int_DDM[:, 0])]
    # Stress_int_point=Stress_int_DDM
    for i in range(Stress_int_DDM.shape[0]):
        Stress_int_point[int(Stress_int_DDM[i,0])-1,:]=Stress_int_DDM[i,:]
    Stress_int_point=np.delete(Stress_int_point, 0, axis=1)

    ## To merge the surface stress output from multiple cpus(Cuong)
    nucleation_data=[]
    for i in range(num_cpus):
        file_elc='./'+foldername_ABAQUS+'/Element_center'+f'{i+1:04d}'+'.txt'
        Surface_elemcenter_temp=np.loadtxt(file_elc)

        file_surf='./'+foldername_ABAQUS+'/Surface_stress'+f'{i+1:04d}'+'.txt'
        if os.path.getsize(file_surf) == 0:
            print("No surface stress in this file.")
        Stress_surface_temp=np.loadtxt(file_surf)
        if Stress_surface_temp.size == 0:
            print("No surface stress in this file.")
        Stress_surface_temp = Stress_surface_temp[:, 1:7]

        nucleation_temp=np.hstack((Surface_elemcenter_temp,Stress_surface_temp))
        if double_write:   # The output from ABAQUS is doubly written in this cycle
            #size_nuclea=len(nucleation_temp)/2
            #(25/07/01 seongjin) slicing number should be integer type
            size_nuclea = int(len(nucleation_temp)/2)
            
            #nucleation_temp=nucleation_temp[:size_nuclea]
            #(25/07/01 seongjin) 2D slicing, not 1D
            nucleation_temp = nucleation_temp[:size_nuclea, :] 
            
        # Deleting contents
        fileID=open(file_elc,'w')   
        fileID.close()     
        fileID=open(file_surf,'w')   
        fileID.close() 
        if len(nucleation_data)==0:
            nucleation_data=nucleation_temp
        else:
            #nucleation_data.extend(nucleation_temp)
            #(25/07/01, seongjin) np.ndarray doesn't have 'extend' 
            nucleation_data = np.vstack([nucleation_data, nucleation_temp])

    nucleation_data=np.array(nucleation_data)
    nucleation_data=nucleation_data[np.argsort(nucleation_data[:, 0])]


    if systemtype_Linux==1:  # For WSL system
        filepath_nu=foldername_temp_wsl+'/'+Jobname+'/'+filename_stress_nucleation
    else:  # For Linux (and WSL2) or MacOS
        filepath_nu='./'+foldername_FEMstress+'/'+filename_stress_nucleation
    fileID = open(filepath_nu, 'w')
    for row in nucleation_data:
        fileID.write(
            f"{int(row[0]):06d} {row[1]:.8e} {row[2]:.8e} {row[3]:.8e} "
            f"{row[4] * FEM_Stress_Scale:.8e} {row[5] * FEM_Stress_Scale:.8e} "
            f"{row[6] * FEM_Stress_Scale:.8e} {row[7] * FEM_Stress_Scale:.8e} "
            f"{row[8] * FEM_Stress_Scale:.8e} {row[9] * FEM_Stress_Scale:.8e}\n"
        )
    fileID.close()

#    print('+++ Finish import stress components.')

    for i in range(num_cpus):
        filename ='./'+foldername_ABAQUS+'/'+filename_FEMoutput_nodes+f'{i+1:04d}'+'.txt'
        fileID = open(filename, 'w')
        fileID.close()

    ## Import link element ID
    link_element=[]
    for i in range(num_cpus):       # Import links found in FEM domain
        file_link='./'+foldername_ABAQUS+'/'+filename_link_element+f'{i+1:04d}'+'.txt'
        link_element_temp=np.loadtxt(file_link)
        fileID = open(file_link, 'w')    # Deleting contents
        fileID.close()
        if i==0:
            link_element=link_element_temp
        else:
            #link_element.extend(link_element_temp)
            #(25/07/01, seongjin) np.ndarray doesn't have 'extend' 
            link_element = np.vstack([link_element, link_element_temp])

    link_element=np.array(link_element)
    for rows in link_element:
        if rows[0]>rows[1]:
            rows[0],rows[1]=rows[1],rows[0]
    link_element=link_element[np.argsort(link_element[:, 0])]

    ## Import rn check
    rn_check=[]
    for i in range(num_cpus):
        file_check='./'+foldername_ABAQUS+'/'+filename_rn_check+f'{i+1:04d}'+'.txt'
        rn_check_temp=np.loadtxt(file_check)
        fileID = open(file_check, 'w')    # Deleting contents
        fileID.close()        
        if i==0:
            rn_check=rn_check_temp
        else:
            #rn_check.extend(rn_check_temp)
            #(25/07/01, seongjin) np.ndarray doesn't have 'extend' 
            rn_check = np.vstack([rn_check, rn_check_temp])

    rn_check=np.array(rn_check)
    rn_check=rn_check[np.argsort(rn_check[:,0])]
    rn_check_length=rn_check.shape[0]

    ## Import element centers
    for i in range(num_cpus):
        file_ele_c='./'+foldername_ABAQUS+'/'+filename_center_element+f'{i+1:04d}'+'.txt'
        fileID = open(file_ele_c, 'w')    # Deleting contents
        fileID.close() 
    
    ## Import rn and links
    rn,links=generate_rnlinks2('./'+foldername_restart+'/restart.data')
    rn,links=np.array(rn),np.array(links)
    #
    length_rn=rn.shape[0]

    for i in range(length_rn):
        if rn[i,3]==0:  # flag-0 node
            j=0
            while j<rn_check_length and rn[i,3]==0:
                if rn[i,4]==rn_check[j,0]:  # Find a match from ABAQUS output
                    if rn_check[j,1]==0:
                        rn[i,3]=2  # this node will be change to flag-0 later
                    elif rn_check[j,1]==7:  # This node is inside impurities
                        rn[i,3]=7  # Fix this node
                    rn_check=np.delete(rn_check,j,axis=0)
                    rn_check_length-=1
                j+=1

        elif rn[i,3]==6:  # For flag-6 node, move the node to its projection of free surface.
            j=0
            while j<rn_check_length and rn[i,3]==6:
                if rn[i,4]==rn_check[j,0]:   # Find a match from ABAQUS output
                    rn[i,3]=2  # this node will be change to flag-0 later.
                    rn_check=np.delete(rn_check,j,axis=0)
                    rn_check_length-=1
                j+=1

        elif rn[i,3]==7:  # For flag-7 node, if it is not found in FEM domain, then change to flag-6
            j=0
            node_inside=False
            while j<rn_check_length and not node_inside:
                if rn[i,4]==rn_check[j,0]: # Find a match from ABAQUS output
                    node_inside=True
                    rn_check=np.delete(rn_check,j,axis=0)
                    rn_check_length-=1
                j+=1
            if not node_inside:  # node not found in FEM domain
                rn[i,3]=6

    for i in range(length_rn):
        if rn[i,3]==0:
            rn[i,3]=6  # this node is outside FEM domain
        elif rn[i,3]==2:
            rn[i,3]=0  # found in FEM domain

    length_links=links.shape[0]

    if length_rn>1:
        # Node numbers in segment(:,1:2)  update (to be the same as the row numbers in rn)
        new_rn_numbers=np.zeros(int(np.max(rn[:,4]))+1,dtype=int)  # this variable stores new node ID 
            # (row number-1 is old ID) new ID=new_rn_numbers(oldID+1)
        for i in range(length_rn):
            new_rn_numbers[int(rn[i,4])]=i+1
        # Node numbers in links update
        for i in range(links.shape[0]):
            links[i,0]=new_rn_numbers[int(links[i,0])]
            links[i,1]=new_rn_numbers[int(links[i,1])]

        for i in range(link_element.shape[0]):
            link_element[i,0]=new_rn_numbers[int(link_element[i,0])]
            link_element[i,1]=new_rn_numbers[int(link_element[i,1])]

        length_link_element=link_element.shape[0]
        # fixing links variable (deleting orphan node)
        i=0
        while i<length_links:  # rm links with just noe node (links from restrat file)
            if links[i,1]==0:
                links=np.delete(links,i,axis=0)
                length_links-=1
                print('    + Broken link detected and deleted (no rn information)')
            else:
                i+=1

        links_Center=np.zeros((length_links,3))
        for i in range(length_links):
            links_Center[i,0]=0.5*(rn[int(links[i,0])-1,0]+rn[int(links[i,1])-1,0])
            links_Center[i,1]=0.5*(rn[int(links[i,0])-1,1]+rn[int(links[i,1])-1,1])
            links_Center[i,2]=0.5*(rn[int(links[i,0])-1,2]+rn[int(links[i,1])-1,2])

#    print(' +++ Done preprocessing ')
    #
    ## Find the element which contains the segment centroid and update segment stress
    #  Update flags of rn if the nodal point is outside the FEM domain.
    if length_rn>1:
        links_element_num=np.zeros((length_links,1)) 
        
        if 'C3D4' in Element_type:  # for tetrahedron elements with four nodes (C3D4 and C3D4T)
            i=0
            while i<length_links and length_link_element>0:
                j=0
                while links_element_num[i]==0 and j<length_link_element:
                    if links[i,0]==link_element[j,0]:
                        if links[i,1]==link_element[j,1]:
                            links_element_num[i]=link_element[j,2]
                            link_element=np.delete(link_element,j,axis=0)
                            length_link_element-=1
                    j+=1
                i+=1
            print(' +++ Done matching segments to elements. ')
        
        elif 'C3D8' in Element_type:  # for hexagonal elements with eight nodes (C3D8, C3D8R, C3D8T, and C3D8RT)
            i=0
            while i<length_links and length_link_element>0:
                j=0
                while links_element_num[i]==0 and j<length_link_element:
                    if links[i,0]==link_element[j,0]:
                        if links[i,1]==link_element[j,1]:
                            links_element_num[i]=link_element[j,2]
                            link_element=np.delete(link_element,j,axis=0)
                            length_link_element-=1
                    j+=1
                i+=1
            # print(' +++ Done matching segments to elements. ')   
    # print(' +++ Done checking nodes.')

## Re-arrange rn, links, and segment, then remove rn, links that are outside of the domain
    if n_step>=1 and length_rn>1:
        # based on flag information (rn(i,4))==6) of two nodal points in
        # links, delete links and then delete rn with 0 neighbor arm #.        
        length_links=links.shape[0]
        if n_step>2:  # Find the liks with on surface nodes
            for i in range(length_links): 
                if rn[int(links[i,0])-1,3]==6 or rn[int(links[i,1])-1,3]==6:
                    node_id=-1
                    if rn[int(links[i,0])-1,3]==6 and rn[int(links[i,1])-1,3]!=6 and rn[int(links[i,1])-1,5]==2:
                        node_id=links[i,1]
                    elif rn[int(links[i,1])-1,3]==6 and rn[int(links[i,0])-1,3]!=6 and rn[int(links[i,0])-1,5]==2:
                        node_id=links[i,0]
                    if node_id>=0:
                        j=0
                        k=0
                        while j<length_links and k<1: 
                            if j!=i:
                                if links[j,0]==node_id and rn[int(links[j,1])-1,3]==6:
                                    rn[int(links[j,0])-1,3]=6
                                    k=k+1
                                elif links[j,1]==node_id and rn[int(links[j,0])-1,3]==6:
                                    rn[int(links[j,1])-1,3]=6
                                    k=k+1
                            j=j+1
        i=0
        while i<length_links: # Remove links with both nodes outside of FEM domain
            if rn[int(links[i,0])-1,3]==6 and rn[int(links[i,1])-1,3]==6:
                rn[int(links[i,0])-1,5]=rn[int(links[i,0])-1,5]-1
                rn[int(links[i,1])-1,5]=rn[int(links[i,1])-1,5]-1
                print('    + Segment '+str(i)+' (Node '+str(links[i,0]-1)+' and Node '+str(links[i,1]-1)+') is outside of the domain: deleted.')
                links=np.delete(links,i,axis=0)
                links_element_num=np.delete(links_element_num,i)
                links_Center=np.delete(links_Center,i,axis=0)
                length_links=length_links-1
            else:
                i=i+1
        length_links=links.shape[0]
        i=0
        while i<length_links: # Remove links debris on the surface
            if links_element_num[i]==0: # link centroid not found in FEM domain
                if rn[int(links[i,0])-1,5]<=1 and rn[int(links[i,1])-1,5]<=1:
                    rn[int(links[i,0])-1,5]=rn[int(links[i,0])-1,5]-1
                    rn[int(links[i,1])-1,5]=rn[int(links[i,1])-1,5]-1
                    print('    + Segment '+str(i)+' (Node '+str(links[i,0]-1)+' and Node '+str(links[i,1]-1)+') is debris segment on surface: deleted.')
                    links=np.delete(links,i,axis=0)
                    links_element_num=np.delete(links_element_num,i)
                    links_Center=np.delete(links_Center,i,axis=0)
                    length_links=length_links-1
                else:
                    i=i+1
            else:
                i=i+1        
        # If a surface node has 1 neighnor node, and the arm's centroid is
        # out side FEM domain, then move the node to the centroid position
        length_links=links.shape[0]
        for i in range(length_links): 
            if links_element_num[i]==0:
                if rn[int(links[i,0])-1,3]==6 and rn[int(links[i,0])-1,5]==1:
                    rn[int(links[i,0])-1,0]=(rn[int(links[i,0])-1,0]+rn[int(links[i,1])-1,0])*0.5
                    rn[int(links[i,0])-1,1]=(rn[int(links[i,0])-1,1]+rn[int(links[i,1])-1,1])*0.5
                    rn[int(links[i,0])-1,2]=(rn[int(links[i,0])-1,2]+rn[int(links[i,1])-1,2])*0.5
                elif rn[int(links[i,1])-1,3]==6 and rn[int(links[i,1])-1,5]==1:
                    rn[int(links[i,1])-1,0]=(rn[int(links[i,0])-1,0]+rn[int(links[i,1])-1,0])*0.5
                    rn[int(links[i,1])-1,1]=(rn[int(links[i,0])-1,1]+rn[int(links[i,1])-1,1])*0.5
                    rn[int(links[i,1])-1,2]=(rn[int(links[i,0])-1,2]+rn[int(links[i,1])-1,2])*0.5

        # If a node has more than 6 neighbors, pin the node (Cuong 1/28/2020)
        for i in range(rn.shape[0]): 
            if rn[i,5]>6: # Node has more than 6 neighbors
                rn[i,3]=7
                k=0
                j=0
                while j<length_links and k<rn[i,5]:
                    if links[j,0]==i+1:
                        rn[int(links[j,1])-1,3]=7
                        k=k+1
                    elif links[j,1]==i+1:
                        rn[int(links[j,0])-1,3]=7
                        k=k+1
                    j=j+1
                print('  Node '+str(i)+' has '+str(int(rn[i,5]))+' neighbors: the node and its neighbors pinned. ')

        length_rn=rn.shape[0]
        for i in range(length_rn): 
            rn[i,4]=i+1 # old id save

        i=0
        while i<length_rn:
            if rn[i,5]<=0:
                rn=np.delete(rn,i,axis=0)
                length_rn=length_rn-1
            else:
                i=i+1

        # sorting rn, links
        new_rn_numbers=np.zeros(int(np.max(rn[:,4])),dtype=int) # this variable stores new node ID
        # (row number-1 is old ID) new ID=new_rn_numbers(oldID+1)
        for i in range(rn.shape[0]):
            new_rn_numbers[int(rn[i,4])-1]=i+1

        for i in range(links.shape[0]):
            links[i,0]=new_rn_numbers[int(links[i,0])-1]
            links[i,1]=new_rn_numbers[int(links[i,1])-1]
        # update links_Center
        #

        ## Write restart.data file
        if systemtype_Linux==1:  # For WSL system
            write_node_data4(foldername_temp_wsl+'/'+Jobname+'/restart',rn,links,xyzlimit)
        else:  # For Linux system or MacOS
            write_node_data4('./'+foldername_restart+'/restart',rn,links,xyzlimit)
        #
    else:
        print(f'    + No nodes exist at Cycle {n_step:d}.')

    ## Write the segment stress (FEMSTRESSXXXXXXXX)
    if systemtype_Linux==1:  # For WSL system
        filepath_stress=foldername_temp_wsl+'/'+Jobname+'/'+filename_FEMstress_head
    else:  # For Linux (and WSL2) or MacOS
        filepath_stress='./'+foldername_FEMstress+'/'+filename_FEMstress_head
    fileID=open(filepath_stress,'w')
    z=0
    uni_stress=500    
    if length_rn>1:
        for i in range(length_links):
            if links_element_num[i]>0:
                if enabled_temp_disp:
                    line_stress='{:06d} {:06d} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.4e}\n'.format(
                        int(links[i,0])-1,int(links[i,1])-1, # Node numbers are adjusted to start from 0 instead of 1.
                        Stress_int_point[int(links_element_num[i])-1,0]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,1]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,2]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,3]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,4]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,5]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,6]
                    )
                else:
                    line_stress='{:06d} {:06d} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e}\n'.format(
                        int(links[i,0])-1,int(links[i,1])-1,
                        Stress_int_point[int(links_element_num[i])-1,0]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,1]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,2]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,3]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,4]*FEM_Stress_Scale,
                        Stress_int_point[int(links_element_num[i])-1,5]*FEM_Stress_Scale
                    )
                fileID.write(line_stress)
            else:
                z=z+1
    else:
        line_stress='{:06d} {:06d} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} 0\n'.format(
            0,0,0.0,0.0,0.0,0.0,0.0,0.0
        )
        fileID.write(line_stress)

    # print(f'   + + + Total {z} segments not assigned.')
    fileID.close()
    ##
    if systemtype_Linux==1:  # For WSL system
        source_dir=os.path.join(foldername_temp_wsl,Jobname)
        dest_dir=os.path.join(foldername_root,foldername_tests,Jobname+'_results',foldername_FEMstress)

        # mv Stress_head file
        src=os.path.join(source_dir,filename_FEMstress_head)
        dst=os.path.join(dest_dir,filename_FEMstress_head)
        shutil.move(src,dst)

        if n_step>=1:
            # mv stress_nucleation file
            src=os.path.join(source_dir,filename_stress_nucleation)
            dst=os.path.join(dest_dir,filename_stress_nucleation)
            shutil.move(src,dst)

            # cp restart and change name
            restart_src=os.path.join(foldername_temp_wsl,Jobname,'restart.data')
            restart_dst_dir=os.path.join(foldername_root,foldername_tests,Jobname+'_results',foldername_restart)
            restart_dst=os.path.join(restart_dst_dir,'restart.data')
            restart_copy=os.path.join(restart_dst_dir,'{:05d}.data'.format(n_step))

            shutil.move(restart_src,restart_dst)
            shutil.copyfile(restart_dst,restart_copy)

# For the segment stress output
# ---------------------------------------------------------------------------------------------
# column  1 : Node ID
# column  2 : Neighbor ID
# column  3 : S11
# column  4 : S22
# column  5 : S33
# column  6 : S23
# column  7 : S13
# column  8 : S12
# column  9 : Flag to determine the segment outside the volume (0: outside,1: inside)
# ---------------------------------------------------------------------------------------------
# For the nucleation site output
# column   1 : RSS
# column   2 : Intg Point coord X
# column   3 : Intg Point coord Y
# column   4 : Intg Point coord Z
# column   5 : Burgers vector X
# column   6 : Burgers vector Y
# column   7 : Burgers vector Z
# column   8 : Slip normal X
# column   9 : Slip normal Y
# column 10 : Slip normal Z
