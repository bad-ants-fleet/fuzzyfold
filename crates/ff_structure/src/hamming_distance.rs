/// Compute Hamming distance between two equal-length strings.
/// This function compares valid UTF-8 strings, therefore it is also case sensitive.
/// Attention: General fucntion, that is not specifically designed for RNA alphabet.
pub fn hamming_distance(seq1: &str, seq2: &str) -> usize {
    assert_eq!(seq1.len(), seq2.len(), "Sequences must have equal length");
    seq1.chars().zip(seq2.chars()).filter(|(x, y)| x != y).count()
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "Sequences must have equal length")]
    fn test_hamming_distance_equal_length() {
        let a: &'static str = "ACUG";
        let b: &'static str = "ACUUU";
        hamming_distance(a, b);
    }

    #[test]
    fn test_hamming_distance_count(){
        let a: &'static str = "ACUG";
        let b: &'static str = "ACUU";
        assert_eq!(hamming_distance(a, b),1);

        // case sensitive
        let lowercase: &'static str = "abc";
        let uppercase: &'static str = "ABC";
        assert_eq!(hamming_distance(lowercase, uppercase), 3);


        // try different alphabet, not just assuming RNA letters
        let c: &'static str = "ACUGYUHU";
        let d: &'static str = "xxUGYUH1";
        assert_eq!(hamming_distance(c, d),3);


        // try hamming distance of two RNA secondary structures in dot-bracket notation
        let a2: &'static str = ".((.....))..";
        let b2: &'static str = "..(((...))).";
        assert_eq!(hamming_distance(a2, b2),5) 



    }
}
