/// Multiformat HMM Parser Demonstration
///
/// This example demonstrates how to use the multiformat HMM parser
/// to parse different HMM formats (HMMER3, PFAM, HMMSearch, InterPro).
///
/// The parser automatically detects the format and parses accordingly.

use omicsx::alignment::MultiFormatHmmParser;

fn main() {
    println!("=== OMICSX Multiformat HMM Parser Demo ===\n");

    // Example 1: HMMER3 Format
    println!("1. Parsing HMMER3 format:");
    let hmmer3_content = r#"HMMER3/f [3.3 | Nov 2019]
NAME  PF00001
ACC   PF00001.28
DESC  7 transmembrane receptor (7tm) superfamily
LENG  345
ALPH  amino
GA    30.00 45.00
TC    35.00 50.00
"#;

    let parser = MultiFormatHmmParser::new();
    match parser.parse_string(hmmer3_content) {
        Ok(profile) => {
            println!("  ✓ Detected: HMMER3");
            println!("  Name: {}", profile.name);
            println!("  Accession: {}", profile.accession.as_deref().unwrap_or("N/A"));
            println!("  Length: {}", profile.length);
            println!("  Alphabet: {}", profile.alphabet);
            println!("  GA Threshold: {}\n", 
                profile.meta.ga_threshold.map_or("None".to_string(), |v| v.to_string()));
        }
        Err(e) => println!("  Error: {}\n", e),
    }

    // Example 2: PFAM Format
    println!("2. Parsing PFAM format:");
    let pfam_content = r#"# STOCKHOLM 1.0
#=GF ID   PF00002
#=GF AC   PF00002.15
#=GF DE   DAG peroxidase superfamily
#=GF LEN  180
#=GF SQ   4589
"#;

    match parser.parse_string(pfam_content) {
        Ok(profile) => {
            println!("  ✓ Detected: PFAM");
            println!("  Name: {}", profile.name);
            println!("  Accession: {}", profile.accession.as_deref().unwrap_or("N/A"));
            println!("  Length: {}\n", profile.length);
        }
        Err(e) => println!("  Error: {}\n", e),
    }

    // Example 3: HMMSearch Format
    println!("3. Parsing HMMSearch format:");
    let hmmsearch_content = r#"# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.3 (Nov 2019); http://hmmer.org/

Query:       PF00003  [M=112]
Accession:   PF00003.12
Description: Zinc finger, C4 type
Scores for complete sequences (score includes all domains):
--- full sequence E-value  score  bias  Description
---   --------- --  -----  -----  ----  -----------
"#;

    match parser.parse_string(hmmsearch_content) {
        Ok(profile) => {
            println!("  ✓ Detected: HMMSearch");
            println!("  Name: {}", profile.name);
            println!("  Accession: {}", profile.accession.as_deref().unwrap_or("N/A"));
            println!("  Length: {}\n", profile.length);
        }
        Err(e) => println!("  Error: {}\n", e),
    }

    // Example 4: InterPro Format
    println!("4. Parsing InterPro format:");
    let interpro_content = r#"ID   IPR000004
AC   IPR000004
DE   Zinc finger, C4-type
DT   28-MAY-2001 (Rel. 7, Created)
DT   19-JAN-2005 (Rel. 16, Last updated)
CC   Short names: ZnF C4; Zn-f_C4
"#;

    match parser.parse_string(interpro_content) {
        Ok(profile) => {
            println!("  ✓ Detected: InterPro");
            println!("  Name: {}", profile.name);
            println!("  Accession: {}", profile.accession.as_deref().unwrap_or("N/A"));
            println!("  Description: {}\n", profile.description.as_deref().unwrap_or("N/A"));
        }
        Err(e) => println!("  Error: {}\n", e),
    }

    // Example 5: List Supported Formats
    println!("5. Supported formats:");
    for format in parser.supported_formats() {
        println!("  • {}", format);
    }
    println!();

    // Example 6: Format Not Recognized
    println!("6. Testing invalid format detection:");
    let invalid_content = "This is not any recognized HMM format";
    match parser.parse_string(invalid_content) {
        Ok(_) => println!("  ✗ Unexpectedly parsed invalid content"),
        Err(e) => println!("  ✓ Correctly rejected: {}\n", e),
    }

    println!("=== Demo Complete ===");
}
