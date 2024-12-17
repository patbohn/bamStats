use std::fs::{create_dir_all, File};
use std::io::{BufWriter, Error, Write};
use std::path::PathBuf;
use std::str::from_utf8;
//use std::io::{BufRead, BufReader};
//use flate2::read::MultiGzDecoder;

use flate2::write::GzEncoder;
use flate2::Compression;

use structopt::StructOpt;

use bam::record::tags::TagValue;
use bam::record::Record;
use bam::BamReader;

///get phred (intergers) to error probability table
pub mod phred_int_to_prob;
use phred_int_to_prob::PHRED_TO_ERROR_PROB;

#[derive(Debug, StructOpt)]
#[structopt(
    name = "bam2fq",
    about = "Calculates QC metrics and extracts tags of unmapped bam output of the dorado basecaller. 
    Also outputs poly A tail estimate located in pt tag and channel and parsed move table."
)]
struct Config {
    #[structopt(parse(from_os_str), help = "Input bam file")]
    input_file: PathBuf,

    #[structopt(
        short = "s",
        long,
        help = "Output directory for stats",
        default_value = "stats_output"
    )]
    output_stats_prefix: String,

    #[structopt(
        short = "f",
        long,
        help = "skip first X nucleotides for quality calculation",
        default_value = "0"
    )]
    skip_first: usize,
}

const MAX_QSCORE: usize = 94;
const MAX_MOVE: usize = 255;
const THREADS: u16 = 8;

fn main() -> std::io::Result<()> {
    let config = Config::from_args();
    get_bam_stats(config)
}

fn get_bam_stats(config: Config) -> Result<(), Error> {
    //define output paths
    let qscore_hist_filename = format!("{}_phred_hist.csv", config.output_stats_prefix);
    let move_hist_filename = format!("{}_move_hist.csv", config.output_stats_prefix);
    let per_read_filename = format!("{}_per_read_stats.csv.gz", config.output_stats_prefix);

    let skip_first: usize = config.skip_first;

    //create output files
    let per_read_statsfile = File::create(per_read_filename)?;
    let mut stats_writer = GzEncoder::new(per_read_statsfile, Compression::fast());
    writeln!(
        &mut stats_writer,
        "read_id_bam,read_id_pod5,read_length,duration,mean_phred,mean_error_rate,poly_a_estimate,mux,channel,read_number,scaling_version,scaling_midpoint,scaling_dispersion,sum_moves,stride"
    )?;

    //read in bam file
    let bam_reader = BamReader::from_path(config.input_file, THREADS).unwrap();
    //read bam let reader = BufReader::new(MultiGzDecoder::new(input_file));

    //create array to count phred qscores and move distribution
    let mut qscore_hist: [u64; MAX_QSCORE] = [0; MAX_QSCORE];
    let mut move_hist: [u64; MAX_MOVE] = [0; MAX_MOVE];

    let mut stride: isize = -1;

    for record in bam_reader {
        //extract read_id, sequence, quality scores per nt, simplex/duplex state and run_id
        let this_record = record.unwrap();

        let bam_id = from_utf8(this_record.name()).unwrap_or("NA");
        let raw_seq = this_record.sequence().to_vec();
        let sequence = String::from_utf8_lossy(&raw_seq);

        let raw_run_id =
            get_optional_string_bam_tag(&this_record, b"RG").unwrap_or("NA".to_owned());
        let run_id = raw_run_id.split("_").next().unwrap_or("NA");
        let poly_a_estimate = get_optional_int_bam_tag(&this_record, b"pt").unwrap_or(-1);
        let mux = get_optional_int_bam_tag(&this_record, b"mx").unwrap_or(-1);
        let channel = get_optional_int_bam_tag(&this_record, b"ch").unwrap_or(-1);
        let read_number = get_optional_int_bam_tag(&this_record, b"rn").unwrap_or(-1);
        let scaling_version =
            get_optional_string_bam_tag(&this_record, b"sv").unwrap_or("NA".to_owned());
        let scaling_midpoint = get_optional_float_bam_tag(&this_record, b"sm").unwrap_or(-1.0);
        let scaling_deviation = get_optional_float_bam_tag(&this_record, b"sd").unwrap_or(-1.0);
        let duration = get_optional_float_bam_tag(&this_record, b"du").unwrap_or(-1.0);
        let parent_id = get_optional_string_bam_tag(&this_record, b"pi").unwrap_or("NA".to_owned());

        let qualities = this_record.qualities().raw();
        for qual in qualities.iter() {
            if qual < &(MAX_QSCORE as u8) {
                qscore_hist[*qual as usize] += 1
            } else {
                panic!("QScore outside of expected range!")
            }

            let moves;

            (stride, moves) = get_nt_durations(&this_record).unwrap();
            for move_duration in moves.iter() {
                if move_duration < &(MAX_MOVE) {
                    move_hist[*move_duration] += 1
                } else {
                    println!("Move duration longer than max: {}", move_duration);
                    move_hist[MAX_MOVE-1] += 1
                }

                let sum_moves = moves.iter().sum::<usize>();

                let mut mean_error_prob = 1.0;
                let mut read_length = 0;
                //(re-)calculate per read mean accuracy (to be seen whether we want to re-evaluate it after trimming?
                if skip_first < qualities.len() {
                    (mean_error_prob, read_length) =
                        calc_mean_median_error(&qualities[..qualities.len() - skip_first]);
                }

                let mean_quality = error_prob_to_phred(mean_error_prob);

                //write per read stats to file
                writeln!(
                    &mut stats_writer,
                    "{},{},{},{:.2},{:.1},{:1.2e},{},{},{},{},{},{},{},{},{}",
                    bam_id,
                    parent_id,
                    read_length,
                    duration,
                    mean_quality,
                    mean_error_prob,
                    poly_a_estimate,
                    mux,
                    channel,
                    read_number,
                    scaling_version,
                    scaling_midpoint,
                    scaling_deviation,
                    sum_moves,
                    stride,
                )?;
            }
        }
    }

    let qual_hist_output_file = File::create(qscore_hist_filename)?;
    let mut writer = BufWriter::new(qual_hist_output_file);
    writeln!(&mut writer, "phred_score,count")?;

    for qscore in 0..MAX_QSCORE {
        writeln!(&mut writer, "{},{}", qscore, qscore_hist[qscore],)?;
    }

    let move_hist_output_file = File::create(move_hist_filename)?;
    let mut writer = BufWriter::new(move_hist_output_file);
    writeln!(&mut writer, "move_duration@stride{},count", stride)?;
    for move_duration in 0..MAX_MOVE {
        writeln!(
            &mut writer,
            "{},{}",
            move_duration, move_hist[move_duration],
        )?;
    }
    Ok(())
}

fn get_optional_string_bam_tag(record: &Record, tag: &[u8; 2]) -> Option<String> {
    if let Some(tag_value) = record.tags().get(tag) {
        match tag_value {
            TagValue::String(tag_string, _) => {
                Some(String::from_utf8_lossy(tag_string).to_string())
            }
            _ => panic!("Unexpected tag type"),
        }
    } else {
        None
    }
}

fn get_optional_int_bam_tag(record: &Record, tag: &[u8; 2]) -> Option<i32> {
    if let Some(tag_value) = record.tags().get(tag) {
        match tag_value {
            TagValue::Int(tag_int, _) => Some(tag_int as i32),
            _ => panic!("Unexpected tag type"),
        }
    } else {
        None
    }
}

fn get_optional_float_bam_tag(record: &Record, tag: &[u8; 2]) -> Option<f32> {
    if let Some(tag_value) = record.tags().get(tag) {
        match tag_value {
            TagValue::Float(tag_float) => Some(tag_float as f32),
            _ => panic!("Unexpected tag type"),
        }
    } else {
        None
    }
}

fn calc_mean_median_error(quality_array: &[u8]) -> (f64, i64) {
    let mut total_prob = 0.0;
    let mut count = 0;

    for quality_score in quality_array.iter() {
        let phred = *quality_score as usize;
        let prob = PHRED_TO_ERROR_PROB[phred];
        total_prob += prob;
        count += 1;
    }
    return (total_prob / count as f64, count as i64);
}

fn error_prob_to_phred(prob: f64) -> f64 {
    return -10.0_f64 * prob.log10();
}

fn phred_to_utf8(quality_array: &[u8]) -> String {
    let mut offset_array: Vec<u8> = Vec::with_capacity(quality_array.len());

    for &byte in quality_array {
        offset_array.push(byte + 33);
    }

    let quality_string = match String::from_utf8(offset_array) {
        Ok(s) => s,
        Err(_e) => {
            panic!("Error converting quality scores to String")
        }
    };
    return quality_string;
}

fn moves_to_utf8(moves: &[usize]) -> (String, Vec<usize>) {
    let mut offset_array: Vec<u8> = Vec::with_capacity(moves.len());
    let mut overflows: Vec<usize> = vec![];

    for &byte in moves {
        if byte < 92 {
            offset_array.push(byte as u8 + 33);
        } else {
            offset_array.push(92 + 33);
            overflows.push(byte);
        }
    }

    let move_string: String = match String::from_utf8(offset_array) {
        Ok(s) => s,
        Err(_e) => {
            panic!("Error converting move durations to String")
        }
    };

    return (move_string, overflows);
}

fn get_nt_durations(record: &Record) -> Option<(isize, Vec<usize>)> {
    if let Some(tag_value) = record.tags().get(b"mv") {
        match tag_value {
            TagValue::IntArray(tag_array) => {
                //println!("Detected mv table");
                Some(calc_moves_per_nt(tag_array.raw()))
            }
            _ => {
                panic!("Unexpected tag type at mv");
            }
        }
    } else {
        None
    }
}

fn calc_moves_per_nt(moves_table: &[u8]) -> (isize, Vec<usize>) {
    //println!("calculating moves per nt");
    let stride = moves_table[0] as isize;

    let mut stay_vector: Vec<usize> = Vec::new();
    let mut count: usize = 1;

    for mv in &moves_table[1..] {
        match mv {
            0 => {
                count += 1;
            }
            1 => {
                stay_vector.push(count);
                count = 1;
            }
            _ => {
                panic!("Unexpected value in move table! Expected only values of <0,1>")
            }
        }
    }
    (stride, stay_vector)
}
