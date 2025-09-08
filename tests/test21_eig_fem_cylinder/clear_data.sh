rm -rf output armdata atomeye dislocation_segs FEMstress properties restart force paradis_output timestep timestep_error
rm -rf *.pyc *.rpy *.rpy.* latest_restart
cp ./ABAQUS/*.inp ./
cp ./ABAQUS/*.for ./
cp ./ABAQUS/*.f ./
rm -rf ./ABAQUS
mkdir ABAQUS
mv ./*.inp ./ABAQUS/  
mv ./*.for ./ABAQUS/  
mv ./*.f ./ABAQUS/  


#rm -rf ./ABAQUS/Data_segs.txt
#rm -rf ./ABAQUS/Data_segs_num.txt
#rm -rf ./ABAQUS/Data_tensor.txt
#rm -rf ./ABAQUS/Restart_Job-0001.odb
#rm -rf ./ABAQUS/VUMAT_DDEM_Nano_20181218.for
#rm -rf ./ABAQUS/abaqus.rpy
#rm -rf ./ABAQUS/res_inp_common1.inp
#rm -rf ./ABAQUS/res_inp_common2.inp
#rm -rf *.ctrl
#rm -rf *.data