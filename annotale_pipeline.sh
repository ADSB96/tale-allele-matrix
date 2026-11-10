#Keep all genomes in one folder. Here, that folder is called "Genomes"
#Changed parentheses to '_'

#step1 (Go to the folder containing the Genomes folder) 

for i in ../Genomes/*.fasta
do
NAM=$(echo $i|cut -f 3 -d '/' | cut -f 1 -d '.')
echo $NAM
java -jar AnnoTALEcli-1.5.jar predict g=$i outdir=$NAM
done

#step2

for i in *
do
if [ -d $i ]
then
java -jar AnnoTALEcli-1.5.jar analyze t=$i/TALE_DNA_sequences.fasta outdir=$i
fi
done

#step3

java -jar AnnoTALEcli-1.5.jar loadAndView outdir=CLASS

#step4

for i in *
do
if [ -f $i/TALE_RVDs*.fasta ]
then
NAM=$(echo $i)
echo $NAM
java -jar AnnoTALEcli-1.5.jar assign c=CLASS/Class_builder_download.xml t=$i/TALE_DNA_parts.fasta s=$STR outdir=$i
fi
done

#step5

for i in */TALE_names*.tsv
do
NAM=$(echo $i|cut -f 1 -d '/')
java -jar AnnoTALEcli-1.5.jar rename r=$i i=$NAM/Genbank__TALE_predictions.gb outdir=$NAM
java -jar AnnoTALEcli-1.5.jar rename r=$i i=$NAM/GFF__TALE_predictions.gff3 outdir=$NAM
done

for i in */TALE_names*.tsv
do
NAM=$(echo $i|cut -f 1 -d '/')
java -jar AnnoTALEcli-1.5.jar rename r=$i i=$NAM/TALE_RVDs.fasta outdir=$NAM
done


############# Create file with all sequences

for i in */Renamed_TALE_sequences/Renamed_TALE_DNA_sequences.fasta
do
NAM=$(echo $i|cut -f 1 -d '/')
sed "s/>/>$NAM|/g" $i >TAL_DNA_$NAM.fa
done

for i in */Renamed_TALE_sequences/Renamed_TALE_protein_sequences.fasta
do
NAM=$(echo $i|cut -f 1 -d '/')
sed "s/>/>$NAM|/g" $i >TAL_AA_$NAM.fa
done

cat TAL_AA* |grep -v '^$' >AFXoc_TAL_AA.fa
cat TAL_DNA* |grep -v '^$' >AFXoc_TAL_DNA.fa

#Added BAI3 as outgroup

#create file with all RVDs

for i in */TALE_RVDs.fasta
do
NAM=$(echo $i|cut -f 1 -d '/')
sed "s/>/>$NAM|/g" $i >TAL_RVD_$NAM.fa
done

cat TAL_RVD* |grep -v '^$' >AFXoc_TAL_RVD.fa


for i in */TALE_names.tsv
do
NAM=$(echo $i|cut -f 1 -d '/')
sed "s/	Tal/	$NAM|Tal/g" $i |sed "s/temp/$NAM|temp/g" | grep -v '#' >>AFXoc_TAL_Names.tsv
done


#Remove repeated strain
#MB-P13-1-1_1
# Removed from files that start with "AF"


## get locations

rm AFXoc_TAL_Locations.tsv
for i in */Renamed_GFF__TALE_loci_1.gff3
do
NAM=$(echo $i|cut -f 1 -d '/')
sed "s/^/$NAM\t/g" $i >>AFXoc_TAL_Locations.tsv
done
