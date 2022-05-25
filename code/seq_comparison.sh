conda activate meme_dep
export PATH=$PATH://home/steveped/meme_install/share:/home/steveped/meme_install/lib:/home/steveped/meme_install/bin

meme-chip -version

# Merge the JASPAR & HOCOMO databases
cp data/external/JASPAR2022_CORE.meme data/external/MOTIFS.meme
sed -n '10,7359p' data/external/HOCOMOCO.v11.meme >> data/external/MOTIFS.meme

# Compare all3 with TFAP2B against all3 without TFAP2B
meme-chip \
  -oc output/meme/all4_tfap2b_vs_notfap2b \
  -db data/external/MOTIFS.meme \
  -neg output/ar_foxa1_gata3_no_tfap2b.fa \
  -desc 'Comparison of binding sites for all TFs based Presence/Absence of TFAP2B' \
  output/ar_foxa1_gata3_tfap2b.fa

# Compare those marked as active regulatory regions against those which are inactive
meme-chip \
  -oc output/meme/all4_h3k27ac_vs_noh3k27ac \
  -db data/external/MOTIFS.meme \
  -neg output/all4_no_h3k27ac.fa \
  -desc 'Comparison of binding sites for all TFs based on overlap with H3K27ac-derived features' \
  output/all4_h3k27ac.fa

# Just analyse the H3K27ac associated ones without performing the comparison
meme-chip \
  -oc output/meme/all4_h3k27ac \
  -db data/external/MOTIFS.meme \
  -desc 'General profile of binding sites for all TFs based on overlap with H3K27ac-derived features' \
  output/all4_h3k27ac.fa

# Try just running FIMO using the core motifs
fimo \
 --oc output/meme/all4_fimo \
 --thresh 1e-2 \
 data/external/all4.meme \
 output/all4_h3k27ac.fa
