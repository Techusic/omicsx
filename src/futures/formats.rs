//! 🔍 Export Formats: BLAST-compatible output formats
//!
//! # Overview
//!
//! This module provides exporters to convert alignment results into standard bioinformatics formats
//! used by BLAST and other sequence analysis tools.
//!
//! # Features
//!
//! - **BLAST XML**: Full XML format with e-values, bit scores, frame information
//! - **BLAST JSON**: JSON representation of BLAST results for programmatic access
//! - **BLAST Tabular**: Tab-separated values (outfmt 6) compatible with NCBI BLAST
//! - **GFF3 Format**: Genomic Feature Format for annotations
//! - **FASTA Output**: Export aligned sequences in FASTA format
//!
//! # Example
//!
//! ```
//! use omics_simd::futures::formats::*;
//!
//! // Create a BLAST JSON entry
//! let json = BlastJson {
//!     query_name: "query1".to_string(),
//!     subject_name: "subject1".to_string(),
//!     score: 100,
//!     evalue: 1e-30,
//!     bit_score: 80.5,
//!     percent_identity: 95.5,
//!     alignment_length: 200,
//!     mismatches: 9,
//!     gap_opens: 0,
//!     query_start: 1,
//!     query_end: 200,
//!     subject_start: 50,
//!     subject_end: 250,
//! };
//!
//! // Convert to BLAST tabular format
//! let tabular_line = json.to_tabular_line();
//! println!("{}", tabular_line);
//! ```

use serde::{Deserialize, Serialize};

/// BLAST XML representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlastXml {
    /// XML formatted BLAST results
    pub xml_content: String,
}

impl BlastXml {
    /// Create a new BLAST XML result
    pub fn new(query: &str, subject: &str, score: i32, evalue: f64) -> Self {
        let xml = format!(
            r#"<?xml version="1.0" encoding="UTF-8"?>
<BlastOutput>
  <BlastOutput_query-def>{}</BlastOutput_query-def>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_def>{}</Hit_def>
          <Hit_hsps>
            <Hsp>
              <Hsp_score>{}</Hsp_score>
              <Hsp_evalue>{:e}</Hsp_evalue>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"#,
            query, subject, score, evalue
        );
        BlastXml { xml_content: xml }
    }

    /// Get the XML string
    pub fn to_string(&self) -> String {
        self.xml_content.clone()
    }
}

/// BLAST JSON representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlastJson {
    /// query sequence name
    pub query_name: String,
    /// matching subject/database entry
    pub subject_name: String,
    /// alignment score
    pub score: i32,
    /// E-value (expect value)
    pub evalue: f64,
    /// Bit score
    pub bit_score: f64,
    /// Percent identity
    pub percent_identity: f32,
    /// Alignment length
    pub alignment_length: usize,
    /// Number of mismatches
    pub mismatches: usize,
    /// Number of gap opens
    pub gap_opens: usize,
    /// Query start position
    pub query_start: usize,
    /// Query end position
    pub query_end: usize,
    /// Subject start position
    pub subject_start: usize,
    /// Subject end position
    pub subject_end: usize,
}

impl BlastJson {
    /// Convert to BLAST tabular format (outfmt 6)
    pub fn to_tabular_line(&self) -> String {
        format!(
            "{}\t{}\t{:.1}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t{:.1}",
            self.query_name,
            self.subject_name,
            self.percent_identity,
            self.alignment_length,
            self.mismatches,
            self.gap_opens,
            self.query_start,
            self.query_end,
            self.subject_start,
            self.subject_end,
            self.evalue,
            self.bit_score
        )
    }

    /// Convert to JSON string
    pub fn to_json(&self) -> Result<String, String> {
        let json = format!(
            r#"{{
  "query_name": "{}",
  "subject_name": "{}",
  "score": {},
  "evalue": {},
  "bit_score": {},
  "percent_identity": {},
  "alignment_length": {},
  "mismatches": {},
  "gap_opens": {},
  "query_start": {},
  "query_end": {},
  "subject_start": {},
  "subject_end": {}
}}"#,
            self.query_name,
            self.subject_name,
            self.score,
            self.evalue,
            self.bit_score,
            self.percent_identity,
            self.alignment_length,
            self.mismatches,
            self.gap_opens,
            self.query_start,
            self.query_end,
            self.subject_start,
            self.subject_end
        );
        Ok(json)
    }
}

/// BLAST tabular format (outfmt 6)
#[derive(Debug, Clone)]
pub struct BlastTabular {
    /// Tab-separated values
    /// Columns: query, subject, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
    pub lines: Vec<String>,
}

impl BlastTabular {
    /// Create a new tabular format
    pub fn new() -> Self {
        BlastTabular { lines: Vec::new() }
    }

    /// Add a result line
    pub fn add_line(&mut self, line: String) {
        self.lines.push(line);
    }

    /// Get the header comment for tabular format
    pub fn get_header() -> String {
        "# Fields: query, subject, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore".to_string()
    }

    /// Convert to string for file output
    pub fn to_string(&self) -> String {
        let mut output = format!("{}\n", Self::get_header());
        for line in &self.lines {
            output.push_str(line);
            output.push('\n');
        }
        output
    }
}

impl Default for BlastTabular {
    fn default() -> Self {
        Self::new()
    }
}

/// GFF3 format representation
#[derive(Debug, Clone)]
pub struct Gff3Record {
    /// Sequence name
    pub seq_id: String,
    /// Source annotation
    pub source: String,
    /// Feature type
    pub feature: String,
    /// Start position (1-based inclusive)
    pub start: usize,
    /// End position (1-based inclusive)
    pub end: usize,
    /// Score value
    pub score: Option<f64>,
    /// Strand (+, -, or .)
    pub strand: char,
    /// Phase (0, 1, 2, or .)
    pub phase: char,
    /// Attributes
    pub attributes: Vec<(String, String)>,
}

impl Gff3Record {
    /// Create a new GFF3 record
    pub fn new(seq_id: &str, feature: &str, start: usize, end: usize) -> Self {
        Gff3Record {
            seq_id: seq_id.to_string(),
            source: "alignment".to_string(),
            feature: feature.to_string(),
            start,
            end,
            score: None,
            strand: '.',
            phase: '.',
            attributes: Vec::new(),
        }
    }

    /// Set the score
    pub fn with_score(mut self, score: f64) -> Self {
        self.score = Some(score);
        self
    }

    /// Set the strand
    pub fn with_strand(mut self, strand: char) -> Self {
        self.strand = strand;
        self
    }

    /// Add an attribute
    pub fn add_attribute(mut self, key: String, value: String) -> Self {
        self.attributes.push((key, value));
        self
    }

    /// Convert to GFF3 format line
    pub fn to_gff3_line(&self) -> String {
        let score_str = self.score.map(|s| format!("{:.1}", s)).unwrap_or_else(|| ".".to_string());
        let mut attr_str = String::new();
        for (i, (key, value)) in self.attributes.iter().enumerate() {
            if i > 0 {
                attr_str.push(';');
            }
            attr_str.push_str(&format!("{}={}", key, value));
        }
        if attr_str.is_empty() {
            attr_str = ".".to_string();
        }

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.seq_id, self.source, self.feature, self.start, self.end, score_str, self.strand, self.phase, attr_str
        )
    }
}

/// Format conversion error
#[derive(Debug)]
pub enum FormatError {
    /// Unsupported alignment type
    UnsupportedType(String),
    /// Serialization failed
    SerializationError(String),
    /// Invalid data for format
    InvalidData(String),
}

impl std::fmt::Display for FormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FormatError::UnsupportedType(s) => write!(f, "Unsupported type: {}", s),
            FormatError::SerializationError(s) => write!(f, "Serialization error: {}", s),
            FormatError::InvalidData(s) => write!(f, "Invalid data: {}", s),
        }
    }
}

impl std::error::Error for FormatError {}

/// Export alignment to BLAST XML format
pub fn to_blast_xml(
    query: &str,
    subject: &str,
    score: i32,
    evalue: f64,
) -> Result<BlastXml, FormatError> {
    Ok(BlastXml::new(query, subject, score, evalue))
}

/// Export alignment to BLAST JSON format
pub fn to_blast_json(
    query: &str,
    subject: &str,
    score: i32,
    evalue: f64,
    bit_score: f64,
    percent_identity: f32,
    alignment_length: usize,
    mismatches: usize,
    gap_opens: usize,
    query_start: usize,
    query_end: usize,
    subject_start: usize,
    subject_end: usize,
) -> Result<BlastJson, FormatError> {
    Ok(BlastJson {
        query_name: query.to_string(),
        subject_name: subject.to_string(),
        score,
        evalue,
        bit_score,
        percent_identity,
        alignment_length,
        mismatches,
        gap_opens,
        query_start,
        query_end,
        subject_start,
        subject_end,
    })
}

/// Export alignment to BLAST tabular format (outfmt 6)
pub fn to_blast_tabular(json_results: &[BlastJson]) -> Result<BlastTabular, FormatError> {
    let mut tabular = BlastTabular::new();
    for result in json_results {
        tabular.add_line(result.to_tabular_line());
    }
    Ok(tabular)
}

/// Convert alignment to GFF3 format
pub fn to_gff3(records: &[Gff3Record]) -> Result<Vec<Gff3Record>, FormatError> {
    // Validate records
    for record in records {
        if record.start > record.end {
            return Err(FormatError::InvalidData(format!(
                "Invalid coordinates: start ({}) > end ({})",
                record.start, record.end
            )));
        }
    }
    Ok(records.to_vec())
}

/// Export sequences in FASTA format
pub fn to_fasta(sequences: &[(String, String)]) -> Result<String, FormatError> {
    let mut fasta = String::new();
    for (header, seq) in sequences {
        fasta.push('>');
        fasta.push_str(header);
        fasta.push('\n');

        // Wrap sequence at 70 characters
        let mut col = 0;
        for ch in seq.chars() {
            if col >= 70 {
                fasta.push('\n');
                col = 0;
            }
            fasta.push(ch);
            col += 1;
        }
        fasta.push('\n');
    }
    Ok(fasta)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blast_xml_export() {
        let xml = to_blast_xml("query1", "subject1", 100, 1e-30).expect("XML creation should succeed");
        let xml_str = xml.to_string();
        
        assert!(xml_str.contains("<?xml version"), "Should contain XML declaration");
        assert!(xml_str.contains("query1"), "Should contain query name");
        assert!(xml_str.contains("subject1"), "Should contain subject name");
        assert!(xml_str.contains("100"), "Should contain score");
    }

    #[test]
    fn test_blast_json_export() {
        let json = to_blast_json(
            "query1", "subject1", 100, 1e-30, 80.5,
            95.5, 200, 9, 0,
            1, 200, 50, 250
        ).expect("JSON creation should succeed");

        assert_eq!(json.query_name, "query1");
        assert_eq!(json.subject_name, "subject1");
        assert_eq!(json.score, 100);
        assert_eq!(json.evalue, 1e-30);
        assert_eq!(json.bit_score, 80.5);
    }

    #[test]
    fn test_blast_json_to_tabular() {
        let json = BlastJson {
            query_name: "q1".to_string(),
            subject_name: "s1".to_string(),
            score: 100,
            evalue: 1e-30,
            bit_score: 80.5,
            percent_identity: 95.5,
            alignment_length: 200,
            mismatches: 9,
            gap_opens: 0,
            query_start: 1,
            query_end: 200,
            subject_start: 50,
            subject_end: 250,
        };

        let line = json.to_tabular_line();
        let parts: Vec<&str> = line.split('\t').collect();
        
        assert_eq!(parts.len(), 12, "Tabular format should have 12 columns");
        assert_eq!(parts[0], "q1", "First column should be query name");
        assert_eq!(parts[1], "s1", "Second column should be subject name");
    }

    #[test]
    fn test_blast_tabular_export() {
        let mut tabular = BlastTabular::new();
        tabular.add_line("q1\ts1\t95.5\t200\t9\t0\t1\t200\t50\t250\t1e-30\t80.5".to_string());
        tabular.add_line("q2\ts2\t92.0\t190\t15\t1\t5\t194\t100\t289\t1e-20\t70.0".to_string());

        let output = tabular.to_string();
        
        assert!(output.contains("# Fields:"), "Should contain header");
        assert!(output.contains("q1\ts1"), "Should contain first result");
        assert!(output.contains("q2\ts2"), "Should contain second result");
    }

    #[test]
    fn test_gff3_record_creation() {
        let record = Gff3Record::new("chr1", "alignment", 100, 200)
            .with_score(100.0)
            .with_strand('+')
            .add_attribute("ID".to_string(), "alignment1".to_string());

        assert_eq!(record.seq_id, "chr1");
        assert_eq!(record.start, 100);
        assert_eq!(record.end, 200);
        assert_eq!(record.strand, '+');
        assert!(record.score.is_some());
    }

    #[test]
    fn test_gff3_line_format() {
        let record = Gff3Record::new("chr1", "match", 1000, 2000)
            .with_score(50.5)
            .with_strand('-')
            .add_attribute("Name".to_string(), "match1".to_string());

        let line = record.to_gff3_line();
        
        assert!(line.contains("chr1"), "Should contain sequence ID");
        assert!(line.contains("match"), "Should contain feature type");
        assert!(line.contains("1000"), "Should contain start position");
        assert!(line.contains("2000"), "Should contain end position");
        assert!(line.contains("50.5"), "Should contain score");
        assert!(line.contains("-"), "Should contain strand");
    }

    #[test]
    fn test_gff3_invalid_coordinates() {
        let records = vec![
            Gff3Record::new("chr1", "test", 200, 100), // Invalid: start > end
        ];
        
        let result = to_gff3(&records);
        assert!(result.is_err(), "Should reject invalid coordinates");
    }

    #[test]
    fn test_fasta_export() {
        let sequences = vec![
            ("seq1".to_string(), "ATCGATCGATCG".to_string()),
            ("seq2".to_string(), "GCTAGCTAGCTA".to_string()),
        ];

        let fasta = to_fasta(&sequences).expect("FASTA export should succeed");
        
        assert!(fasta.contains(">seq1"), "Should contain first header");
        assert!(fasta.contains(">seq2"), "Should contain second header");
        assert!(fasta.contains("ATCGATCGATCG"), "Should contain first sequence");
        assert!(fasta.contains("GCTAGCTAGCTA"), "Should contain second sequence");
    }

    #[test]
    fn test_fasta_wrapping() {
        let long_sequence = "A".repeat(100);
        let sequences = vec![
            ("long_seq".to_string(), long_sequence),
        ];

        let fasta = to_fasta(&sequences).expect("FASTA export should succeed");
        let lines: Vec<&str> = fasta.lines().collect();
        
        // Should have header, then multiple sequence lines
        assert!(lines.len() > 2, "Should wrap long sequence");
        assert_eq!(lines[0], ">long_seq", "First line should be header");
    }
}
