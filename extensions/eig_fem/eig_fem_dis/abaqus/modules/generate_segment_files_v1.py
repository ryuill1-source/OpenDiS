import os
import shutil 

def generate_segment_files_v1(foldername_ABAQUS,foldername_segs,
                    filename_seg_head,filename_seg_temp,filename_seg_lastCycle,
                    filename_seg_num,dt_FEM_each_step,EffSlip,InelasticStep,i):
    # Working folder: ./tests/'Job name'
    # Copy Data_segs.txt files to ABAQUS folder
    # filename_segment=[foldername_segs,'/',filename_seg_head,num2str(i,'%04d')];
    # filename_gradient=[foldername_segs,'/',filename_Grad_head,num2str(i,'%04d')];
    filename_segment=foldername_segs+'/'+filename_seg_head+'_all'
    N_segments=0
    if InelasticStep or i>1:
        if os.path.isfile(filename_segment) :
            shutil.copyfile(filename_segment,foldername_ABAQUS+'/'+filename_seg_temp)
            with open(foldername_ABAQUS+'/'+filename_seg_temp,'r') as f_1:
                for line in f_1:
                    line=line.strip()
                    if line:
                        N_segments+=1
    # Segment information from last ParaDiS cycle
    filename_segment=foldername_segs+'/'+filename_seg_head+'.final'
    N_segments_lastCycle=0
    if InelasticStep or i>1:
        if os.path.isfile(filename_segment) :
            shutil.copyfile(filename_segment,foldername_ABAQUS+'/'+filename_seg_lastCycle)
            with open(foldername_ABAQUS+'/'+filename_seg_lastCycle,'r') as f_2:
                for line in f_2:
                    line=line.strip()
                    if line:
                        N_segments_lastCycle+=1
    
    # Generate Data_segs_num.txt file to provide the number of segments in the
    # current step
    filename=foldername_ABAQUS+'/'+filename_seg_num
    fID=open(filename,'w')
    temp_line = f"{int(N_segments)} {int(N_segments_lastCycle)} {float(dt_FEM_each_step):.5e} {float(EffSlip):.5e}\n"
    fID.write(temp_line)
    fID.close()
    return N_segments