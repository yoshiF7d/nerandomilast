mkdir log/velop
mkdir result
#sbatch -o log/velop/d0.log -J velop_d0 -p Medium --mem 64G --wrap "python3 velopepi.py --name d0 --samples d0"
#sbatch -o log/velop/d8.log -J velop_d8 -p Medium --mem 64G --wrap "python3 velopepi.py --name d8 --samples d8_1 d8_2"
sbatch -o log/velop/d10_BI.log -J velop_d10_BI -p Medium --mem 64G --wrap "python3 velopepi.py --name d10_BI --samples d10_BI_1 d10_BI_2"
sleep 1
sbatch -o log/velop/d10_Cont.log -J velop_d10_Cont -p Medium --mem 64G --wrap "python3 velopepi.py --name d10_Cont --samples d10_Cont_1 d10_Cont_2"
sleep 1
sbatch -o log/velop/d14_BI.log -J velop_d14_BI -p Medium --mem 64G --wrap "python3 velopepi.py --name d14_BI --samples d14_BI_1 d14_BI_2"
sleep 1
sbatch -o log/velop/d14_Cont.log -J velop_d14_Cont -p Medium --mem 64G --wrap "python3 velopepi.py --name d14_Cont --samples d14_Cont_1 d14_Cont_2"
sleep 1
sbatch -o log/velop/d21_BI.log -J velop_d21_BI -p Medium --mem 64G --wrap "python3 velopepi.py --name d21_BI --samples d21_BI_1 d21_BI_2"
sleep 1
sbatch -o log/velop/d21_Cont.log -J velop_d21_Cont -p Medium --mem 64G --wrap "python3 velopepi.py --name d21_Cont --samples d21_Cont_1 d21_Cont_2"
