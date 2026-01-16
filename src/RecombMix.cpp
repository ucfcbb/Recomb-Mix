//  * --------------------------------------------------------------------------------------------------------
//  * Name: RecombMix.cpp
//  * Description: Local Ancestry Inference based on improved Loter with new model, using recombination rate.
//  * Author: Yuan Wei 
//  * Created on: Jan 21, 2023
//  * Modified on: Dec 22, 2025
//  * --------------------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <map>
#include <random>
#include <iomanip>
#include <vector>
#include <queue>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cmath>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <omp.h>
using namespace std;

//remove spaces from string
string remove_spaces(string input_string){
    input_string.erase(remove_if (input_string.begin(), input_string.end(), [](char c){ return (c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == '\v' || c == '\f'); }), input_string.end());
    return input_string;
}

//split a string into tokens with given delimiter
vector<string> split_string(string line, char delimiter){
    vector<string> tokens;
    string token;
    stringstream line_stream(line);
    while (getline(line_stream, token, delimiter)){
        if (!token.empty()){
            tokens.push_back(token);
        }
    }
    return tokens;
}

static void show_usage(string program_name){
    cout << "Usage: " << program_name << " <Option(s)>\n";
    cout << "Option(s):\n";
    cout << "\t-h,--help\t\t\t\t\tShow this help message\n";
    cout << "\t-p,--panel <INPUT PANEL FILE>\t\t\tInput panel path and file name\n";
    cout << "\t-q,--query <INPUT QUERY FILE>\t\t\tInput query path and file name\n";
    cout << "\t-g,--genetic <INPUT GENETIC MAPPING FILE>\tInput genetic mapping path and file name\n";
    cout << "\t-a,--ancestry <INPUT POPULATION ANCESTRY FILE>\tInput population ancestry path and file name\n";
    cout << "\t-o,--output <OUTPUT DIRECTORY PATH>\t\tOutput directory path for all files\n";
    cout << "\t-i,--inferred <OUTPUT INFERRED FILE NAME>\tOutput inferred local ancestry file name\n";
    cout << "\t-e,--weight <WEIGHT>\t\t\t\tWeight of recombination rate in cost function\n";
    cout << "\t-f,--frequency <ALLELE FREQUENCY>\t\tMinor allele frequency threshold\n";
    cout << "\t-u,--outputcompactpanel <IDENTIFIER>\t\tSpecify whether output compact reference panel (1: true; 0: false)\n";
    cout << "\t-s,--estimate <MAXIMUM GAP PHYSICAL DISTANCE>\tMaximum gap physical distance (number of sites) for local ancestry estimation\n";
    cout << "\t-t,--threads <NUMBER OF THREADS>\t\tNumber of threads\n";
}

int main(int argc, char *argv[]){
    try {
        cout << "start program" << "\n";

        //input variables with default values
        string chromosome_id = "."; //chromosome id
        bool has_chromosome_id = false;
        string output_directory_path = "./"; //default is current folder
        double weight = 1.5; //recombination rate weight in cost function
        double maf_threshold = 0.0; //minor allele frequency threshold to include allele values in the graph

        //other variables
        string version_number = "0.7"; //program version number
        int number_of_populations = 0; //number of ways of admixture population
        int number_of_site_values = 2; //biallelic (0 or 1)
        char population_delimiter = '\t'; //text format
        char vcf_delimiter = '\t'; //VCF format
        char map_delimiter = '\t'; //HapMap format

        unordered_map<string, int> sample_populations; //key: sample_id, value: population_id
        unordered_map<int, int> sample_index_populations; //key: sample_index_id, value: population_id
        map<int, double> genetic_maps; //genetic map read from input HapMap data; key: physical position (base pairs), value: genetic position (centiMorgans)
        vector<double> panel_genetic_positions;
        vector<double> panel_recombination_rates;
        vector<double> panel_recombination_rate_penalties; //reciprocal of the normalized recombination rate in [0, 1]; the large recombination rate, the high probability having recombination event, the low penalty
        vector<unordered_map<int, double>> population_zero_frequency_per_site; //population allele frequencies in reference panel; each site has a map having key: population id; value: allele frequency of site having zero value

        vector<bool> query_intersected_sites; //store intersected sites information only
        int number_of_query_individuals = 0;
        int number_of_query_sites = 0;
        vector<string> query_individual_ids;
        vector<int> query_physical_positions;

        int number_of_panel_individuals = 0;
        int number_of_panel_sites = 0;
        vector<string> panel_individual_ids;
        vector<int> panel_physical_positions;

        string input_reference_panel_path_and_file_name;
        string input_query_panel_path_and_file_name;
        string input_map_path_and_file_name;
        string input_population_path_and_file_name;
        string output_inference_individuals_path_and_file_name;
        string output_compact_panel_path_and_file_name;
        string output_compact_panel_population_label_path_and_file_name;
        string output_inference_individuals_file_name = "admix_inferred_ancestral_values_local.txt";
        bool output_compact_panel = false;
        bool has_queries = false;
        int input_buffer_size = 1;
        int maximum_gap_physical_distance = 0;
        unsigned long number_of_threads = 1;

        //get command line arguments
        string program_name = argv[0];
        string program_arguments = "";
        for (int i = 1; i < argc; i++){
            string argument = argv[i];
            if ((argument == "-h") || (argument == "--help")){
                show_usage(argv[0]);
                cout << "end program" << "\n";
                return 0;
            }
            else if ((argument == "-p") || (argument == "--panel")){
                if (i + 1 < argc){
                    input_reference_panel_path_and_file_name = argv[i + 1];
                    program_arguments += " p=" + input_reference_panel_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-p (or --panel) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-q") || (argument == "--query")){
                if (i + 1 < argc){
                    input_query_panel_path_and_file_name = argv[i + 1];
                    has_queries = true;
                    program_arguments += " q=" + input_query_panel_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-q (or --query) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-g") || (argument == "--genetic")){
                if (i + 1 < argc){
                    input_map_path_and_file_name = argv[i + 1];
                    program_arguments += " g=" + input_map_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-g (or --genetic) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-a") || (argument == "--ancestry")){
                if (i + 1 < argc){
                    input_population_path_and_file_name = argv[i + 1];
                    program_arguments += " a=" + input_population_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-a (or --ancestry) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-o") || (argument == "--output")){
                if (i + 1 < argc){
                    output_directory_path = argv[i + 1];
                    program_arguments += " o=" + output_directory_path;
                    i++;
                }
                else {
                    cout << "-o (or --output) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-i") || (argument == "--inferred")){
                if (i + 1 < argc){
                    output_inference_individuals_file_name = argv[i + 1];
                    program_arguments += " i=" + output_inference_individuals_file_name;
                    i++;
                }
                else {
                    cout << "-i (or --inferred) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-e") || (argument == "--weight")){
                if (i + 1 < argc){
                    weight = stod(argv[i + 1]);
                    program_arguments += " e=" + to_string(weight);
                    i++;
                }
                else {
                    cout << "-e (or --weight) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-f") || (argument == "--frequency")){
                if (i + 1 < argc){
                    maf_threshold = stod(argv[i + 1]);
                    program_arguments += " f=" + to_string(maf_threshold);
                    i++;
                }
                else {
                    cout << "-f (or --frequency) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-u") || (argument == "--outputcompactpanel")){
                if (i + 1 < argc){
                    int output_compact_panel_int = stoi(argv[i + 1]);
                    output_compact_panel = output_compact_panel_int == 1 ? true : false;
                    program_arguments += " u=" + to_string(output_compact_panel_int);
                    i++;
                }
                else {
                    cout << "-u (or --outputcompactpanel) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-s") || (argument == "--estimate")){
                if (i + 1 < argc){
                    maximum_gap_physical_distance = stoi(argv[i + 1]);
                    program_arguments += " s=" + to_string(maximum_gap_physical_distance);
                    i++;
                }
                else {
                    cout << "-s (or --estimate) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-t") || (argument == "--threads")){
                if (i + 1 < argc){
                    number_of_threads = stoul(argv[i + 1]);
                    program_arguments += " t=" + to_string(number_of_threads);
                    i++;
                }
                else {
                    cout << "-t (or --threads) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else {
                cout << "unrecognized arguments: " << argument << "\n";
                show_usage(argv[0]);
                cout << "end program" << "\n";
                return 1;
            }
        }
        program_arguments = program_arguments.size() > 0 ? program_arguments.substr(1) : "";

        //construct output files
        if (output_directory_path.substr(output_directory_path.size() - 1).compare("/") == 0){
            output_directory_path = output_directory_path.substr(0, output_directory_path.size() - 1);
        }
        output_inference_individuals_path_and_file_name = output_directory_path + "/" + output_inference_individuals_file_name;
        output_compact_panel_path_and_file_name = output_directory_path + "/compact_reference_panel.vcf.gz";
        output_compact_panel_population_label_path_and_file_name = output_directory_path + "/compact_reference_panel_population_label.txt";

        //output variables
        cout << "program name: " << program_name << "\nprogram arguments: " << program_arguments << "\n";
        cout << "weight=" << fixed << setprecision(6) << weight << "\n";
        unsigned long number_of_max_threads = omp_get_max_threads();
        if (number_of_threads < 1 || number_of_threads > number_of_max_threads){
            number_of_threads = number_of_max_threads;
        }
        cout << "number_of_threads=" << number_of_threads << endl;
        cout << "maximum_gap_physical_distance=" << maximum_gap_physical_distance << endl;
        cout << "output files:" << "\n";
        cout << "local ancestray inference: " << output_inference_individuals_path_and_file_name << "\n";
        if (output_compact_panel){
            cout << "compact reference panel: " << output_compact_panel_path_and_file_name << "\n";
            cout << "compact reference panel population label: " << output_compact_panel_population_label_path_and_file_name << "\n";
        }
        
        //set vcf file header index; vcf file should have: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT (9 fields) before the first individual id
        int chromosome_id_location_on_file = 0;
        int position_location_on_file = 1;
        int id_location_on_file = 2;
        int reference_base_location_on_file = 3;
        int alternate_base_location_on_file = 4;
        int quality_location_on_file = 5;
        int filter_location_on_file = 6;
        int information_location_on_file = 7;
        int format_location_on_file = 8;
        int individual_id_start_position = 9;

        //variables for input and output
        ifstream input_file_data;
        ofstream output_file_data;

        cout << "start read reference ancestry data" << "\n";

        //use two maps as bidirectional map to store origin labels and ids
        unordered_map<string, int> origin_to_ids;
        unordered_map<int, string> origin_to_labels;
        int label_id = 0;
        set<string> population_labels;

        //read sample population label from file
        input_file_data.open(input_population_path_and_file_name);
        if (!input_file_data){
            cout << "cannot open file " + input_population_path_and_file_name << "\n";
            exit(1);
        }
        if (input_file_data.is_open()){
            int line_number = 0;
            string line_from_file;
            while (getline(input_file_data, line_from_file)){
                if (line_from_file.substr(0, 1).compare("#") == 0){
                    //do nothing; skip the header or comment line
                }
                else {
                    vector<string> tokens = split_string(line_from_file, population_delimiter);
                    if (tokens.size() >= 2){
                        //population format: sample_id \t population_label
                        string sample_id = remove_spaces(tokens[0]);
                        string population_label = remove_spaces(tokens[1]);
                        if (origin_to_ids.count(population_label) <= 0){
                            origin_to_ids[population_label] = label_id;
                            origin_to_labels[label_id] = population_label;
                            label_id++;
                        }
                        sample_populations[sample_id] = origin_to_ids[population_label];
                        sample_index_populations[line_number] = origin_to_ids[population_label];
                        population_labels.insert(population_label);
                    }
                    line_number++;
                }
            }
            input_file_data.close();
        }
        number_of_populations = population_labels.size();
        cout << "number_of_populations=" << number_of_populations << "\n";

        cout << "end read reference ancestry data" << "\n";

        cout << "start read panel data" << "\n";

        //read panel data from file
        streamsize io_buffer_size = input_buffer_size * 1024 * 1024; //default is 1 MB buffer size
        ifstream input_panel_file;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> stream_in_buffer_panel;
        if (input_reference_panel_path_and_file_name.substr(input_reference_panel_path_and_file_name.size() - 2, 2).compare("gz") == 0){
            input_panel_file.open(input_reference_panel_path_and_file_name, ios_base::in | ios_base::binary);
            stream_in_buffer_panel.push(boost::iostreams::gzip_decompressor());
        }
        else {
            input_panel_file.open(input_reference_panel_path_and_file_name, ios_base::in);
        }
        if (!input_panel_file){
            cout << "cannot open file " + input_reference_panel_path_and_file_name << "\n";
            exit(1);
        }
        stream_in_buffer_panel.push(input_panel_file, io_buffer_size);
        istream input_vcf_file_data(&stream_in_buffer_panel);
        if (input_panel_file.is_open()){
            int line_number = 0;
            bool end_of_file = false;
            string line_from_file;
            string line_from_file_next;
            while (!end_of_file){
                line_from_file = line_from_file_next;
                if (!end_of_file){
                    getline(input_vcf_file_data, line_from_file_next);
                    if (input_vcf_file_data.eof()){
                        end_of_file = true;
                    }
                }
                if (line_number == 0){
                    line_number++;
                    continue;
                }
                if (line_from_file.substr(0, 2).compare("##") == 0){
                    //do nothing; skip header lines
                }
                else if (line_from_file.substr(0, 1).compare("#") == 0){
                    //get panel individual ids
                    vector<string> tokens = split_string(line_from_file, vcf_delimiter);
                    for (int i = individual_id_start_position; i < tokens.size(); i++){
                        panel_individual_ids.push_back(tokens[i]);
                    }
                    number_of_panel_individuals = tokens.size() - individual_id_start_position;
                }
                else {
                    //population allele frequency calculation (key: population id; value: a pair of count of alleles having zero value and count of total alleles)
                    unordered_map<int, tuple<int, int>> population_zero_count_total;

                    //get panel site values and population allele frequency of the site
                    vector<bool> site_of_panel_individuals;
                    vector<string> tokens = split_string(line_from_file, vcf_delimiter);
                    for (int i = 0; i < tokens.size(); i++){
                        //get physical position
                        if (i == position_location_on_file){
                            panel_physical_positions.push_back(stoi(tokens[i]));
                        }
                        //get chromosome id
                        if (!has_chromosome_id && i == chromosome_id_location_on_file){
                            chromosome_id = tokens[i];
                            has_chromosome_id = true;
                        }
                        //get panel individual first and second haplotype site values and allele frequency of the site
                        if (tokens[i].substr(1, 1).compare("|") == 0 || tokens[i].substr(1, 1).compare("/") == 0){ //phased or unphased
                            int site_of_first_haplotype = stoi(tokens[i].substr(0, 1));
                            int site_of_second_haplotype = stoi(tokens[i].substr(2, 1));
                            site_of_panel_individuals.push_back(site_of_first_haplotype);
                            site_of_panel_individuals.push_back(site_of_second_haplotype);

                            //build compact panel
                            int individual_id = i - individual_id_start_position;
                            int haplotype_1_id = individual_id * 2;
                            int haplotype_2_id = haplotype_1_id + 1;
                            string sample_id = panel_individual_ids[individual_id];
                            int population_label_id = -1;
                            if (sample_populations.size() > 0 && sample_populations.count(sample_id) > 0){
                                population_label_id = sample_populations[sample_id];
                            }
                            else {
                                population_label_id = -1; //"UNKNOWN";
                            }
                            if (population_label_id > -1){
                                if (population_zero_count_total.count(population_label_id) <= 0){
                                    population_zero_count_total[population_label_id] = tuple<int, int>(0, 0);
                                }
                                if (site_of_first_haplotype == 0){
                                    get<0>(population_zero_count_total[population_label_id]) += 1;
                                }
                                get<1>(population_zero_count_total[population_label_id]) += 1;
                                if (site_of_second_haplotype == 0){
                                    get<0>(population_zero_count_total[population_label_id]) += 1;
                                }
                                get<1>(population_zero_count_total[population_label_id]) += 1;
                            }
                        }
                    }
                    if (site_of_panel_individuals.size() > 0){
                        site_of_panel_individuals.clear();

                        //store population allele frequency of the site (key: population id; value: population allele frequency of site whose value is zero (note: population allele frequency of site whose value is one = 1 - population allele frequency of site whose value is zero))
                        unordered_map<int, double> population_zero_frequency;
                        for (const auto & [key_population_label_id, value] : population_zero_count_total){
                            population_zero_frequency[key_population_label_id] = double(get<0>(value)) / double(get<1>(value));
                        }
                        population_zero_frequency_per_site.push_back(population_zero_frequency);

                        number_of_panel_sites++;
                    }
                    else {
                        cout << "data at line " + to_string(line_number) + " is not valid" << "\n";
                        cout << line_from_file << "\n";
                    }
                }
                line_number++;
            }
            boost::iostreams::close(stream_in_buffer_panel);
            input_panel_file.close();
        }
        
        cout << "end read panel data" << "\n";

        //output compact reference panel
        if (output_compact_panel){
            ofstream output_panel_file(output_compact_panel_path_and_file_name, ios_base::out | ios_base::binary);
            if (!output_panel_file){
                cout << "cannot create or open file " + output_compact_panel_path_and_file_name << "\n";
                exit(1);
            }
            boost::iostreams::filtering_streambuf<boost::iostreams::output> stream_out_buffer;
            stream_out_buffer.push(boost::iostreams::gzip_compressor());
            stream_out_buffer.push(output_panel_file, io_buffer_size);
            ostream output_vcf_file_data(&stream_out_buffer);
            if (output_panel_file.is_open()){
                output_vcf_file_data << "##fileformat=VCFv4.2" << "\n";
                output_vcf_file_data << "##source=RecombMix_v" << version_number << "\n";
                output_vcf_file_data << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << "\n";
                output_vcf_file_data << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
                for (int i = 0; i < origin_to_labels.size(); i++){
                    output_vcf_file_data << origin_to_labels[i];
                    if (i < origin_to_labels.size() - 1){
                        output_vcf_file_data << "\t";
                    }
                    else {
                        output_vcf_file_data << "\n";
                    }
                }
            }
            for (int i = 0; i < population_zero_frequency_per_site.size(); i++){
                vector<string> individuals_site_value(origin_to_labels.size());
                for (const auto & [key_population_label_id, value_allele_frequency_zero] : population_zero_frequency_per_site[i]){
                    string site_value_str = "";
                    if (value_allele_frequency_zero > 0.0){
                        if (1.0 - value_allele_frequency_zero > 0.0){
                            site_value_str = "0/1";
                        }
                        else {
                            site_value_str = "0/0";
                        }
                    }
                    else {
                        site_value_str = "1/1";
                    }
                    individuals_site_value[key_population_label_id] = site_value_str;
                }
                if (output_panel_file.is_open()){
                    output_vcf_file_data << chromosome_id << "\t" << panel_physical_positions[i] << "\t.\t.\t.\t.\tPASS\t.\tGT\t";
                    for (int j = 0; j < individuals_site_value.size(); j++){
                        output_vcf_file_data << individuals_site_value[j];
                        if (j < individuals_site_value.size() - 1){
                            output_vcf_file_data << "\t";
                        }
                        else {
                            output_vcf_file_data << "\n";
                        }
                    }
                }
            }
            boost::iostreams::close(stream_out_buffer);
            output_panel_file.close();

            output_file_data.open(output_compact_panel_population_label_path_and_file_name, ios::trunc);
            if (!output_file_data){
                cout << "cannot create or open file " + output_compact_panel_population_label_path_and_file_name << "\n";
                exit(1);
            }
            if (output_file_data.is_open()){
                output_file_data << "#Sample_ID\tPopulation_Label" << "\n";
                for (int i = 0; i < origin_to_labels.size(); i++){
                    output_file_data << origin_to_labels[i] << "\t" << origin_to_labels[i] << "\n";
                }
            }
            output_file_data.close();
        }

        if (!has_queries){
            cout << "end program" << "\n";
            return 0;
        }

        cout << "start read query data" << "\n";

        //read query data from file
        int physical_position_panel_current = 0;
        vector<int> panel_physical_position_index_filtered;
        vector<int> query_physical_position_index_filtered;
        ifstream input_query_file;
        boost::iostreams::filtering_streambuf<boost::iostreams::input> stream_in_buffer_query;
        if (input_query_panel_path_and_file_name.substr(input_query_panel_path_and_file_name.size() - 2, 2).compare("gz") == 0){
            input_query_file.open(input_query_panel_path_and_file_name, ios_base::in | ios_base::binary);
            stream_in_buffer_query.push(boost::iostreams::gzip_decompressor());
        }
        else {
            input_query_file.open(input_query_panel_path_and_file_name, ios_base::in);
        }
        if (!input_query_file){
            cout << "cannot open file " + input_query_panel_path_and_file_name << "\n";
            exit(1);
        }
        stream_in_buffer_query.push(input_query_file, io_buffer_size);
        istream input_query_file_data(&stream_in_buffer_query);
        if (input_query_file.is_open()){
            int line_number = 0;
            bool end_of_file = false;
            string line_from_file;
            string line_from_file_next;
            while (!end_of_file){
                line_from_file = line_from_file_next;
                if (!end_of_file){
                    getline(input_query_file_data, line_from_file_next);
                    if (input_query_file_data.eof()){
                        end_of_file = true;
                    }
                }
                if (line_number == 0){
                    line_number++;
                    continue;
                }
                
                if (line_from_file.substr(0, 2).compare("##") == 0){
                    //do nothing; skip header lines
                }
                else if (line_from_file.substr(0, 1).compare("#") == 0){
                    //get query individual ids
                    vector<string> tokens = split_string(line_from_file, vcf_delimiter);
                    for (int i = individual_id_start_position; i < tokens.size(); i++){
                        query_individual_ids.push_back(tokens[i]);
                    }
                    number_of_query_individuals = tokens.size() - individual_id_start_position;
                }
                else {
                    //get query site values
                    bool is_query_site_intersected = false;
                    vector<bool> site_of_query_individuals;
                    vector<string> tokens = split_string(line_from_file, vcf_delimiter);
                    for (int i = 0; i < tokens.size(); i++){
                        //get physical position
                        if (i == position_location_on_file){
                            int physical_position_query = stoi(tokens[i]);
                            query_physical_positions.push_back(physical_position_query);
                            while (physical_position_panel_current < panel_physical_positions.size() && panel_physical_positions[physical_position_panel_current] < physical_position_query){
                                physical_position_panel_current++;
                            }
                            if (physical_position_panel_current < panel_physical_positions.size() && physical_position_query == panel_physical_positions[physical_position_panel_current]){
                                panel_physical_position_index_filtered.push_back(physical_position_panel_current);
                                query_physical_position_index_filtered.push_back(number_of_query_sites);
                                is_query_site_intersected = true;
                            }
                            else {
                                //this site from query is not in the reference panel, skip it
                                number_of_query_sites++;
                                break;
                            }
                        }
                        //get query individual first and second haplotype site values
                        if (i >= individual_id_start_position){
                            if (tokens[i].substr(1, 1).compare("|") == 0){ //phased biallelic value
                                int site_of_first_haplotype = stoi(tokens[i].substr(0, 1));
                                int site_of_second_haplotype = stoi(tokens[i].substr(2, 1));
                                site_of_query_individuals.push_back(site_of_first_haplotype);
                                site_of_query_individuals.push_back(site_of_second_haplotype);
                            }
                            else {
                                site_of_query_individuals.clear();
                                break;
                            }
                        }
                    }
                    if (site_of_query_individuals.size() > 0 && is_query_site_intersected){
                        //store current site values from query individuals
                        query_intersected_sites.insert(query_intersected_sites.end(), site_of_query_individuals.begin(), site_of_query_individuals.end());
                        site_of_query_individuals.clear();
                        number_of_query_sites++;
                    }
                }
                line_number++;
            }
            boost::iostreams::close(stream_in_buffer_query);
            input_query_file.close();
        }

        cout << "end read query data" << "\n";

        size_t number_of_intersected_sites = panel_physical_position_index_filtered.size();
        cout << "number_of_intersecting_markers=" << number_of_intersected_sites;
        cout << ", percentage_of_markers_in_query_covered_by_reference=" << fixed << setprecision(2) << (query_physical_positions.size() > 0 ? (double(panel_physical_position_index_filtered.size()) / double(query_physical_positions.size())) * 100 : 0) << "%" << "\n";

        cout << "start read map data" << "\n";

        //read genetic map data from file
        input_file_data.open(input_map_path_and_file_name);
        if (!input_file_data){
            cout << "cannot open file " + input_map_path_and_file_name << "\n";
            exit(1);
        }
        if (input_file_data.is_open()){
            int line_number = 0;
            string line_from_file;
            while (getline(input_file_data, line_from_file)){
                //get genetic map data
                vector<string> tokens = split_string(line_from_file, map_delimiter);
                if (tokens.size() >= 4){
                    //genetic map HapMap format: #Chromosome \t Position(bp) \t Rate(cM/Mb) \t Map(cM)
                    if (line_number > 0){
                        //skip the first header line
                        int physical_location = stoi(tokens[1]);
                        double genetic_location = stod(tokens[3]);
                        genetic_maps[physical_location] = genetic_location;
                    }
                }
                line_number++;
            }
            input_file_data.close();
        }

        //interpolate genetic map based on panel physical positions (filtered index)
        map<int, double>::iterator position_exact_iterator, position_lower_iterator, position_upper_iterator;
        for (int i = 0; i < panel_physical_position_index_filtered.size(); i++){
            double panel_genetic_position = 0.0;
            //search current physical position from genetic map
            if (!genetic_maps.empty()){
                position_exact_iterator = genetic_maps.find(panel_physical_positions[panel_physical_position_index_filtered[i]]);
                if (position_exact_iterator == genetic_maps.end()){
                    //exact physical position not found; need to estimate genetic position with neighbor positions in genetic map
                    position_upper_iterator = genetic_maps.upper_bound(panel_physical_positions[panel_physical_position_index_filtered[i]]);
                    if (position_upper_iterator == genetic_maps.end()){
                        if (panel_physical_positions[panel_physical_position_index_filtered[i]] > genetic_maps.rbegin()->first){
                            //physical position is larger than any physical position in map
                            panel_genetic_position = genetic_maps.rbegin()->second;
                        }
                        else { 
                            //physical position is less than any physical position in map
                            panel_genetic_position = genetic_maps.begin()->second;
                        }
                    }
                    else {
                        if (position_upper_iterator == genetic_maps.begin()){
                            //physical position is less than any physical position in map
                            panel_genetic_position = genetic_maps.begin()->second;
                        }
                        else {
                            //estimate genetic position by linear interpolation
                            int physical_position_next = position_upper_iterator->first;
                            double genetic_position_next = position_upper_iterator->second;
                            position_lower_iterator = --position_upper_iterator;
                            int physical_position_previous = position_lower_iterator->first;
                            double genetic_position_previous = position_lower_iterator->second;
                            panel_genetic_position = genetic_position_next + (double)(panel_physical_positions[panel_physical_position_index_filtered[i]] - physical_position_next) * (genetic_position_next - genetic_position_previous) / (double)(physical_position_next - physical_position_previous);
                        }
                    }
                }
                else {
                    //exact physical position found; just add its corrlated genetic position
                    panel_genetic_position = position_exact_iterator->second;
                }
            }
            panel_genetic_positions.push_back(panel_genetic_position);
        }

        //calculate recombination rates (HapMap genetic map)
        double rr_min = INFINITY;
        double rr_max = numeric_limits<double>::min();
        double rr_curr = 0.0;
        for (int i = 0; i < panel_genetic_positions.size(); i++){
            if (i + 1 < panel_genetic_positions.size()){
                //the current position (or i)'s recombination rate is cM/Mb, which is cM/(bp/(10^6)), i.e. r[i]=(g[i+1]-g[i])/((p[i+1]-p[i])/(10^6))
                rr_curr = (panel_genetic_positions[i + 1] - panel_genetic_positions[i]) / (double(panel_physical_positions[panel_physical_position_index_filtered[i + 1]] - panel_physical_positions[panel_physical_position_index_filtered[i]]) / pow(10, 6));
            }
            else {
                //the last position's recombination rate is zero
                rr_curr = 0.0;
            }
            if (rr_curr < rr_min){
                rr_min = rr_curr;
            }
            if (rr_curr > rr_max){
                rr_max = rr_curr;
            }
            panel_recombination_rates.push_back(rr_curr);
        }

        //calculate recombination rate penalties
        double rrp_curr = 0.0;
        for (int i = 0; i < panel_recombination_rates.size(); i++){
            //normalize recombination rate into [0, 1], and then add 1, becoming [1, 2], then take its reciprocal, becoming [0.5, 1.0], then times 2 to make it within the range of [1.0, 2.0]
            rrp_curr = (rr_max - rr_min > 0 && panel_recombination_rates[i] - rr_min >= 0) ? (panel_recombination_rates[i] - rr_min) / (rr_max - rr_min) : 0.0;
            panel_recombination_rate_penalties.push_back(2.0 / (rrp_curr + 1));
        }

        cout << "end read map data" << "\n";

        cout << "start verify data" << "\n";

        bool is_data_satisfied = true;

        //verify parameters
        if (panel_physical_position_index_filtered.size() <= 0){
            is_data_satisfied = false;
            cout << "reference panel has no intersecting markers with query" << "\n";
        }
        if (number_of_populations == 0){
            is_data_satisfied = false;
            cout << "reference panel has no ancestry information" << "\n";
        }
        if (weight < 0.0){
            is_data_satisfied = false;
            cout << "weight of recombination rate is invalid" << "\n";
        }

        cout << "end verify data" << "\n";

        if (!is_data_satisfied){
            cout << "end program" << "\n";
            return 0;
        }

        int number_of_query_haplotypes = number_of_query_individuals * 2;
        int number_of_panel_haplotypes = number_of_panel_individuals * 2;
        cout << "number of query haplotypes: " << to_string(number_of_query_haplotypes) << "\n";
        cout << "number of query sites: " << to_string(number_of_query_sites) << "\n";
        cout << "number of panel haplotypes: " << to_string(number_of_panel_haplotypes) << "\n";
        cout << "number of panel sites: " << to_string(number_of_panel_sites) << "\n";

        cout << "start run queries" << "\n";

        output_file_data.open(output_inference_individuals_path_and_file_name, ios::trunc);
        if (!output_file_data){
            cout << "cannot create or open file " + output_inference_individuals_path_and_file_name << "\n";
            exit(1);
        }
        if (output_file_data.is_open()){
            output_file_data << "#Population label and ID: ";
            for (int i = 0; i < origin_to_labels.size(); i++){
                output_file_data << origin_to_labels[i] << "=" << i;
                if (i < origin_to_labels.size() - 1){
                    output_file_data << "\t";
                }
                else {
                    output_file_data << "\n";
                }
            }
        }
        output_file_data.close();
        vector<string> query_output_buffer(number_of_query_haplotypes);

        #pragma omp parallel for
        for (int query_haplotype_index = 0; query_haplotype_index < number_of_query_haplotypes; query_haplotype_index++){
            vector<double> scores(number_of_intersected_sites * (number_of_populations * 2), INFINITY);
            vector<int> node_ids_with_min_score_per_population(number_of_intersected_sites * number_of_populations, -1);
            vector<int> node_ids_with_min_score_overall(number_of_intersected_sites, -1);
            vector<int> paths(number_of_intersected_sites * (number_of_populations * 2), -1);

            //scan sites in reverse order
            for (int site_id = 0; site_id < number_of_intersected_sites; site_id++){
                for (auto & [key_population_label_id, value_allele_frequency_zero] : population_zero_frequency_per_site[panel_physical_position_index_filtered[number_of_intersected_sites - 1 - site_id]]){
                    //apply minor allele frequency threshold
                    if (value_allele_frequency_zero <= maf_threshold){
                        value_allele_frequency_zero = 0.0;
                    }
                    else if ((1.0 - value_allele_frequency_zero) <= maf_threshold){
                        value_allele_frequency_zero = 1.0;
                    }
                    else {
                        value_allele_frequency_zero = 0.5;
                    }

                    //calculate mismatch penalty score
                    double allele_value_zero_mismatch_penalty = INFINITY;
                    double allele_value_one_mismatch_penalty = INFINITY;
                    int query_site_value = query_intersected_sites[(number_of_intersected_sites - 1 - site_id) * number_of_query_haplotypes + query_haplotype_index];
                    if (value_allele_frequency_zero == 0.5){
                        if (query_site_value == 0){
                            allele_value_zero_mismatch_penalty = 0.0;
                            allele_value_one_mismatch_penalty = 1.0;
                        }
                        else {
                            allele_value_zero_mismatch_penalty = 1.0;
                            allele_value_one_mismatch_penalty = 0.0;
                        }
                    }
                    else if (value_allele_frequency_zero == 0.0){
                        if (query_site_value == 0){
                            allele_value_one_mismatch_penalty = 1.0;
                        }
                        else {
                            allele_value_one_mismatch_penalty = 0.0;
                        }
                    }
                    else {
                        if (query_site_value == 0){
                            allele_value_zero_mismatch_penalty = 0.0;
                        }
                        else {
                            allele_value_zero_mismatch_penalty = 1.0;
                        }
                    }

                    //calculate total penalty score
                    if (site_id == 0){
                        scores[0 * (number_of_populations * 2) + (key_population_label_id * 2)] = allele_value_zero_mismatch_penalty;
                        scores[0 * (number_of_populations * 2) + (key_population_label_id * 2 + 1)] = allele_value_one_mismatch_penalty;
                        if (allele_value_zero_mismatch_penalty > allele_value_one_mismatch_penalty){
                            node_ids_with_min_score_per_population[0 * number_of_populations + key_population_label_id] = key_population_label_id * 2 + 1;
                            if (node_ids_with_min_score_overall[0] == -1 || scores[0 * (number_of_populations * 2) + node_ids_with_min_score_overall[0]] > scores[0 * (number_of_populations * 2) + (key_population_label_id * 2 + 1)]){
                                node_ids_with_min_score_overall[0] = key_population_label_id * 2 + 1;
                            }
                        }
                        else {
                            node_ids_with_min_score_per_population[0 * number_of_populations + key_population_label_id] = key_population_label_id * 2;
                            if (node_ids_with_min_score_overall[0] == -1 || scores[0 * (number_of_populations * 2) + node_ids_with_min_score_overall[0]] > scores[0 * (number_of_populations * 2) + (key_population_label_id * 2)]){
                                node_ids_with_min_score_overall[0] = key_population_label_id * 2;
                            }
                        }
                    }
                    else {
                        //calculate penalty score for both population nodes
                        double template_change_penalty = weight * panel_recombination_rate_penalties[number_of_intersected_sites - 1 - site_id];
                        if (scores[(site_id - 1) * (number_of_populations * 2) + node_ids_with_min_score_per_population[(site_id - 1) * number_of_populations + key_population_label_id]] > scores[(site_id - 1) * (number_of_populations * 2) + node_ids_with_min_score_overall[site_id - 1]] + template_change_penalty){
                            scores[site_id * (number_of_populations * 2) + (key_population_label_id * 2)] = allele_value_zero_mismatch_penalty + scores[(site_id - 1) * (number_of_populations * 2) + node_ids_with_min_score_overall[site_id - 1]] + template_change_penalty;
                            paths[site_id * (number_of_populations * 2) + (key_population_label_id * 2)] = node_ids_with_min_score_overall[site_id - 1];
                            scores[site_id * (number_of_populations * 2) + (key_population_label_id * 2 + 1)] = allele_value_one_mismatch_penalty + scores[(site_id - 1) * (number_of_populations * 2) + node_ids_with_min_score_overall[site_id - 1]] + template_change_penalty;
                            paths[site_id * (number_of_populations * 2) + (key_population_label_id * 2 + 1)] = node_ids_with_min_score_overall[site_id - 1];
                        }
                        else {
                            scores[site_id * (number_of_populations * 2) + (key_population_label_id * 2)] = allele_value_zero_mismatch_penalty + scores[(site_id - 1) * (number_of_populations * 2) + node_ids_with_min_score_per_population[(site_id - 1) * number_of_populations + key_population_label_id]];
                            paths[site_id * (number_of_populations * 2) + (key_population_label_id * 2)] = node_ids_with_min_score_per_population[(site_id - 1) * number_of_populations + key_population_label_id];
                            scores[site_id * (number_of_populations * 2) + (key_population_label_id * 2 + 1)] = allele_value_one_mismatch_penalty + scores[(site_id - 1) * (number_of_populations * 2) + node_ids_with_min_score_per_population[(site_id - 1) * number_of_populations + key_population_label_id]];
                            paths[site_id * (number_of_populations * 2) + (key_population_label_id * 2 + 1)] = node_ids_with_min_score_per_population[(site_id - 1) * number_of_populations + key_population_label_id];
                        }

                        if (node_ids_with_min_score_per_population[site_id * number_of_populations + key_population_label_id] < 0 || scores[site_id * (number_of_populations * 2) + node_ids_with_min_score_per_population[site_id * number_of_populations + key_population_label_id]] > scores[site_id * (number_of_populations * 2) + (key_population_label_id * 2)]){
                            node_ids_with_min_score_per_population[site_id * number_of_populations + key_population_label_id] = key_population_label_id * 2;
                        }
                        if (node_ids_with_min_score_per_population[site_id * number_of_populations + key_population_label_id] < 0 || scores[site_id * (number_of_populations * 2) + node_ids_with_min_score_per_population[site_id * number_of_populations + key_population_label_id]] > scores[site_id * (number_of_populations * 2) + (key_population_label_id * 2 + 1)]){
                            node_ids_with_min_score_per_population[site_id * number_of_populations + key_population_label_id] = key_population_label_id * 2 + 1;
                        }

                        if (node_ids_with_min_score_overall[site_id] < 0 || scores[site_id * (number_of_populations * 2) + node_ids_with_min_score_overall[site_id]] > scores[site_id * (number_of_populations * 2) + (key_population_label_id * 2)]){
                            node_ids_with_min_score_overall[site_id] = key_population_label_id * 2;
                        }
                        if (node_ids_with_min_score_overall[site_id] < 0 || scores[site_id * (number_of_populations * 2) + node_ids_with_min_score_overall[site_id]] > scores[site_id * (number_of_populations * 2) + (key_population_label_id * 2 + 1)]){
                            node_ids_with_min_score_overall[site_id] = key_population_label_id * 2 + 1;
                        }
                    }
                }
            }

            //find node with min score for the first site
            int node_id_selected = 0;
            for (int node_id = 1; node_id < number_of_populations * 2; node_id++){
                if (scores[(number_of_intersected_sites - 1) * (number_of_populations * 2) + node_id_selected] > scores[(number_of_intersected_sites - 1) * (number_of_populations * 2) + node_id]){
                    node_id_selected = node_id;
                }
            }

            //trace back population labels while merge and output consecutive sites having the same ancestral value
            if (maximum_gap_physical_distance <= 0){
                //merge consecutive sites having the same ancestral value
                query_output_buffer[query_haplotype_index] = query_individual_ids[query_haplotype_index / 2] + "_" + to_string(query_haplotype_index % 2);
                int start_position_prev = -2;
                int population_label_prev = -2;
                int filtered_site_id = number_of_intersected_sites - 1;
                for (int site_id = 0; site_id < number_of_query_sites; site_id++){
                    int population_label_id_selected = -1;
                    if (filtered_site_id >= 0 && query_physical_position_index_filtered[number_of_intersected_sites - 1 - filtered_site_id] == site_id){
                        population_label_id_selected = node_id_selected / 2;
                        if (filtered_site_id > 0){
                            node_id_selected = paths[filtered_site_id * (number_of_populations * 2) + node_id_selected];
                            filtered_site_id -= 1;
                        }
                    }
                    if (population_label_id_selected != population_label_prev){
                        if (population_label_prev != -2){
                            query_output_buffer[query_haplotype_index] += "\t" + to_string(start_position_prev) + "\t" + to_string(query_physical_positions[site_id - 1]) + "\t" + to_string(population_label_prev);
                        }
                        start_position_prev = query_physical_positions[site_id];
                        population_label_prev = population_label_id_selected;
                    }
                }
                query_output_buffer[query_haplotype_index] += "\t" + to_string(start_position_prev) + "\t" + to_string(query_physical_positions[number_of_query_sites - 1]) + "\t" + to_string(population_label_prev) + "\n";
            }
            else {
                //merge consecutive sites having the same ancestral value while estimate the gap ancestry label
                query_output_buffer[query_haplotype_index] = query_individual_ids[query_haplotype_index / 2] + "_" + to_string(query_haplotype_index % 2);
                int start_site_id_second_last = -1;
                int population_label_id_second_last = -1;
                int start_site_id_last = -1;
                int population_label_id_last = -1;
                int filtered_site_id = number_of_intersected_sites - 1;
                for (int site_id = 0; site_id < number_of_query_sites; site_id++){
                    int population_label_id_selected = -1;
                    if (filtered_site_id >= 0 && query_physical_position_index_filtered[number_of_intersected_sites - 1 - filtered_site_id] == site_id){
                        population_label_id_selected = node_id_selected / 2;
                        if (filtered_site_id > 0){
                            node_id_selected = paths[filtered_site_id * (number_of_populations * 2) + node_id_selected];
                            filtered_site_id -= 1;
                        }
                    }
                    if (population_label_id_last > -1){
                        if (population_label_id_selected > -1){
                            if (population_label_id_last != population_label_id_selected){
                                query_output_buffer[query_haplotype_index] += "\t" + to_string(query_physical_positions[start_site_id_last]) + "\t" + to_string(query_physical_positions[site_id - 1]) + "\t" + to_string(population_label_id_last);
                                population_label_id_last = population_label_id_selected;
                                start_site_id_last = site_id;
                            }
                        }
                        else {
                            population_label_id_second_last = population_label_id_last;
                            start_site_id_second_last = start_site_id_last;
                            population_label_id_last = -1;
                            start_site_id_last = site_id;
                        }
                    }
                    else {
                        if (population_label_id_second_last > -1){
                            if (population_label_id_selected > -1){
                                if ((population_label_id_second_last == population_label_id_selected) && (site_id - start_site_id_last <= maximum_gap_physical_distance)){
                                    population_label_id_last = population_label_id_second_last;
                                    start_site_id_last = start_site_id_second_last;
                                    population_label_id_second_last = -1;
                                }
                                else {
                                    query_output_buffer[query_haplotype_index] += "\t" + to_string(query_physical_positions[start_site_id_second_last]) + "\t" + to_string(query_physical_positions[start_site_id_last - 1]) + "\t" + to_string(population_label_id_second_last);
                                    population_label_id_second_last = -1;
                                    population_label_id_last = population_label_id_selected;
                                    start_site_id_last = site_id;
                                }
                            }
                        }
                        else {
                            if (population_label_id_selected > -1){
                                population_label_id_last = population_label_id_selected;
                                start_site_id_last = site_id;
                            }
                        }
                    }
                }
                if (population_label_id_second_last > -1){
                    query_output_buffer[query_haplotype_index] += "\t" + to_string(query_physical_positions[start_site_id_second_last]) + "\t" + to_string(query_physical_positions[start_site_id_last - 1]) + "\t" + to_string(population_label_id_second_last);
                }
                if (population_label_id_last > -1){
                    query_output_buffer[query_haplotype_index] += "\t" + to_string(query_physical_positions[start_site_id_last]) + "\t" + to_string(query_physical_positions[number_of_query_sites - 1]) + "\t" + to_string(population_label_id_last);
                }
            }
        }

        cout << "end run queries" << "\n";

        output_file_data.open(output_inference_individuals_path_and_file_name, ios::app);
        if (!output_file_data){
            cout << "cannot create or open file " + output_inference_individuals_path_and_file_name << "\n";
            exit(1);
        }
        if (output_file_data.is_open()){
            //output file format: query_individual_id_query_haplotype_index \t origin_label_ids of each segment's start position, end position, population value (merged if consecutive segments having the same ancestral value) separated by \t
            for (int query_haplotype_index = 0; query_haplotype_index < number_of_query_haplotypes; query_haplotype_index++){   
                output_file_data << query_output_buffer[query_haplotype_index] << endl;
            }
        }
        output_file_data.close();

        cout << "end program" << "\n";
    }
    catch (exception & exception_output){
        cerr << exception_output.what() << "\n";
    }
    return 0;
}