%Chk=G2.chk
%nprocshared=8
%mem=1GB
# HF/6-31G** SCRF=(Solvent=Water) POP=MK iop(6/50=1) geom=allcheck guess=read

G2.esp
