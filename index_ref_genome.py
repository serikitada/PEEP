mkdir tmp
java -Xmx4g -jar FlashFry-assembly-1.12.jar \
 index \
 --tmpLocation ./tmp \
 --database chr22_cas9ngg_database \
 --reference chr22.fa.gz \
 --enzyme spcas9ngg