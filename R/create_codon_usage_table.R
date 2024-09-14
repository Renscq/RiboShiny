#' Function to create the codon table.
#'
#' Codon table contains the following columns:
#' \itemize{} 
#' \item Codon: the codon sequence
#' \item AntiCodon: the anticodon sequence
#' \item AA: the amino acid
#' \item Abbreviation: the abbreviation of the amino acid
#' \item Name: the name of the amino acid
#' 
#' @param seq sequence type, either "DNA" or "RNA"
#'
#' @return Returns a \code{Dataframe} object.
#'
#' @examples
#'
#' codon_table <- create_codon_usage_table(seq = "DNA")
#' codon_table <- create_codon_usage_table(seq = "DNA")
#'
#' @export
# 


create_codon_usage_table <- function(
    seq = "DNA"
    ){
  
  Codon = c('TGA', 'TAA', 'TAG', 'GCT', 'GCA', 'GCC', 'GCG', 'TGT', 'TGC', 'GAT', 'GAC', 'GAA', 'GAG', 'TTT', 'TTC', 'GGA', 
            'GGT', 'GGG', 'GGC', 'CAT', 'CAC', 'ATT', 'ATA', 'ATC', 'AAA', 'AAG', 'TTG', 'CTT', 'TTA', 'CTA', 'CTG', 'CTC', 
            'ATG', 'AAT', 'AAC', 'CCA', 'CCT', 'CCC', 'CCG', 'CAA', 'CAG', 'AGA', 'AGG', 'CGT', 'CGA', 'CGG', 'CGC', 'TCT', 
            'TCA', 'AGT', 'AGC', 'TCC', 'TCG', 'ACT', 'ACA', 'ACC', 'ACG', 'GTT', 'GTG', 'GTA', 'GTC', 'TGG', 'TAT', 'TAC')
  
  AntiCodon = c('TCA', 'TTA', 'CTA', 'ACG', 'TGC', 'GGC', 'CGC', 'ACA', 'GCA', 'ATC', 'GTC', 'TTC', 'CTC', 'AAA', 'GAA', 'TTC', 
                'ACC', 'CCC', 'GCC', 'ATA', 'GTG', 'AAT', 'TAT', 'GAT', 'TTT', 'CTT', 'CAA', 'AAG', 'TAA', 'TAG', 'CAG', 'GAG', 
                'CAT', 'TTA', 'GTT', 'TGG', 'AGG', 'GGG', 'CGG', 'TTG', 'CTG', 'TCT', 'CCT', 'ACG', 'TCG', 'CCG', 'GCG', 'AGA', 
                'AGA', 'ACT', 'AGC', 'GGA', 'CGA', 'TGA', 'TGT', 'GGT', 'CGT', 'AAC', 'CAC', 'TAC', 'GAC', 'CCA', 'ATA', 'GTA')
  
  AA <- c('TER', 'TER', 'TER', 'Ala', 'Ala', 'Ala', 'Ala', 'Cys', 'Cys', 'Asp', 'Asp', 'Glu', 'Glu', 'Phe', 'Phe', 'Gly', 
          'Gly', 'Gly', 'Gly', 'His', 'His', 'Ile', 'Ile', 'Ile', 'Lys', 'Lys', 'Leu', 'Leu', 'Leu', 'Leu', 'Leu', 'Leu', 
          'Met', 'Asn', 'Asn', 'Pro', 'Pro', 'Pro', 'Pro', 'Gln', 'Gln', 'Arg', 'Arg', 'Arg', 'Arg', 'Arg', 'Arg', 'Ser', 
          'Ser', 'Ser', 'Ser', 'Ser', 'Ser', 'Thr', 'Thr', 'Thr', 'Thr', 'Val', 'Val', 'Val', 'Val', 'Trp', 'Tyr', 'Tyr')

  Abbr = c('*', '*', '*', 'A', 'A', 'A', 'A', 'C', 'C', 'D', 'D', 'E', 'E', 'F', 'F', 'G',
           'G', 'G', 'G', 'H', 'H', 'I', 'I', 'I', 'K', 'K', 'L', 'L', 'L', 'L', 'L', 'L',
           'M', 'N', 'N', 'P', 'P', 'P', 'P', 'Q', 'Q', 'R', 'R', 'R', 'R', 'R', 'R', 'S',
           'S', 'S', 'S', 'S', 'S', 'T', 'T', 'T', 'T', 'V', 'V', 'V', 'V', 'W', 'Y', 'Y')
  
  Name = c('Termination', 'Termination', 'Termination', 'Alanine', 'Alanine', 'Alanine', 'Alanine', 'Cysteine', 'Cysteine',
           'Aspartic Acid', 'Aspartic Acid', 'Glutamic Acid', 'Glutamic Acid', 'Phenylalanine', 'Phenylalanine', 'Glycine',
           'Glycine', 'Glycine', 'Glycine', 'Histidine', 'Histidine', 'Isoleucine', 'Isoleucine', 'Isoleucine', 'Lysine', 
           'Lysine', 'Leucine', 'Leucine', 'Leucine', 'Leucine', 'Leucine', 'Leucine', 'Methionine', 'Asparagine', 
           'Asparagine', 'Proline', 'Proline', 'Proline', 'Proline', 'Glutamine', 'Glutamine', 'Arginine', 'Arginine', 
           'Arginine', 'Arginine', 'Arginine', 'Arginine', 'Serine', 'Serine', 'Serine', 'Serine', 'Serine', 'Serine', 
           'Threonine', 'Threonine', 'Threonine', 'Threonine', 'Valine', 'Valine', 'Valine', 'Valine', 'Tryptophan',
           'Tyrosine', 'Tyrosine')
  
  codon_table <- data.frame(
    Codon = Codon,
    AntiCodon = AntiCodon,
    AA = AA,
    Abbr = Abbr,
    Name = Name
  )
  
  # convert the DNA to RNA sequence
  if (seq == "RNA") {
    codon_table$Codon <- gsub("T", "U", codon_table$Codon)
    codon_table$AntiCodon <- gsub("T", "U", codon_table$AntiCodon)
  }

  return(codon_table)
}