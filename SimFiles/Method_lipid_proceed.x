#!/bin/bash
RUNS=$1
for i in `seq 1 $RUNS`
do
 mkdir Run_$i
 cd Run_$i
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod_2.in -o Production_0_100.out -c ../Store/Solv_0_100.inpcrd -p ../Store/Solv_0_100.prmtop -r 0_100.restart -x 0_100.mdcrd
 fail=0
 until grep -q "Total wall time" Production_0_100.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod_2.in -o Production_0_100.out -c ../Store/Solv_0_100.inpcrd -p ../Store/Solv_0_100.prmtop -r 0_100.restart -x 0_100.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_10_90.out -c 0_100.restart -p ../Store/Solv_10_90.prmtop -r 10_90.restart -x 10_90.mdcrd
 fail=0
 until grep -q "Total wall time" Production_10_90.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod_2.in -o Production_0_100.out -c ../Store/Solv_0_100.inpcrd -p ../Store/Solv_0_100.prmtop -r 0_100.restart -x 0_100.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_10_90.out -c 0_100.restart -p ../Store/Solv_10_90.prmtop -r 10_90.restart -x 10_90.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_20_80.out -c 10_90.restart -p ../Store/Solv_20_80.prmtop -r 20_80.restart -x 20_80.mdcrd
 fail=0
 until grep -q "Total wall time" Production_20_80.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_10_90.out -c 0_100.restart -p ../Store/Solv_10_90.prmtop -r 10_90.restart -x 10_90.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_20_80.out -c 10_90.restart -p ../Store/Solv_20_80.prmtop -r 20_80.restart -x 20_80.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_30_70.out -c 20_80.restart -p ../Store/Solv_30_70.prmtop -r 30_70.restart -x 30_70.mdcrd
 fail=0
 until grep -q "Total wall time" Production_30_70.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_20_80.out -c 10_90.restart -p ../Store/Solv_20_80.prmtop -r 20_80.restart -x 20_80.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_30_70.out -c 20_80.restart -p ../Store/Solv_30_70.prmtop -r 30_70.restart -x 30_70.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_40_60.out -c 30_70.restart -p ../Store/Solv_40_60.prmtop -r 40_60.restart -x 40_60.mdcrd
 fail=0
 until grep -q "Total wall time" Production_40_60.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_30_70.out -c 20_80.restart -p ../Store/Solv_30_70.prmtop -r 30_70.restart -x 30_70.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_40_60.out -c 30_70.restart -p ../Store/Solv_40_60.prmtop -r 40_60.restart -x 40_60.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_50_50.out -c 40_60.restart -p ../Store/Solv_50_50.prmtop -r 50_50.restart -x 50_50.mdcrd
 fail=0
 until grep -q "Total wall time" Production_50_50.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_40_60.out -c 30_70.restart -p ../Store/Solv_40_60.prmtop -r 40_60.restart -x 40_60.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_50_50.out -c 40_60.restart -p ../Store/Solv_50_50.prmtop -r 50_50.restart -x 50_50.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_60_40.out -c 50_50.restart -p ../Store/Solv_60_40.prmtop -r 60_40.restart -x 60_40.mdcrd
 fail=0
 until grep -q "Total wall time" Production_60_40.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_50_50.out -c 40_60.restart -p ../Store/Solv_50_50.prmtop -r 50_50.restart -x 50_50.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_60_40.out -c 50_50.restart -p ../Store/Solv_60_40.prmtop -r 60_40.restart -x 60_40.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_70_30.out -c 60_40.restart -p ../Store/Solv_70_30.prmtop -r 70_30.restart -x 70_30.mdcrd
 fail=0
 until grep -q "Total wall time" Production_70_30.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_60_40.out -c 50_50.restart -p ../Store/Solv_60_40.prmtop -r 60_40.restart -x 60_40.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_70_30.out -c 60_40.restart -p ../Store/Solv_70_30.prmtop -r 70_30.restart -x 70_30.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_80_20.out -c 70_30.restart -p ../Store/Solv_80_20.prmtop -r 80_20.restart -x 80_20.mdcrd
 fail=0
 until grep -q "Total wall time" Production_80_20.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_70_30.out -c 60_40.restart -p ../Store/Solv_70_30.prmtop -r 70_30.restart -x 70_30.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_80_20.out -c 70_30.restart -p ../Store/Solv_80_20.prmtop -r 80_20.restart -x 80_20.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_90_10.out -c 80_20.restart -p ../Store/Solv_90_10.prmtop -r 90_10.restart -x 90_10.mdcrd
 fail=0
 until grep -q "Total wall time" Production_90_10.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_80_20.out -c 70_30.restart -p ../Store/Solv_80_20.prmtop -r 80_20.restart -x 80_20.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_90_10.out -c 80_20.restart -p ../Store/Solv_90_10.prmtop -r 90_10.restart -x 90_10.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_100_0.out -c 90_10.restart -p ../Store/Solv_100_0.prmtop -r 100_0.restart -x 100_0.mdcrd
 fail=0
 until grep -q "Total wall time" Production_100_0.out
 do
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_90_10.out -c 80_20.restart -p ../Store/Solv_90_10.prmtop -r 90_10.restart -x 90_10.mdcrd
  nohup pmemd.cuda_SPFP -O -i ~/Python/SimFiles/Mutation/Lipid_Prod.in -o Production_100_0.out -c 90_10.restart -p ../Store/Solv_100_0.prmtop -r 100_0.restart -x 100_0.mdcrd
  fail=$(($fail+1))
  if [ $fail -eq 10 ]
  then
   break
  fi
 done
 cd ../
done
