//  * --------------------------------------------------------------------------------------------------------
//  * Name: generate_compact_panel.cpp
//  * Description: Create a compact reference panel VCF file and an ancestry label file, from an original reference panel VCF file and an ancestry label file.
//  * Author: Yuan Wei 
//  * Created on: Jul 03, 2026
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
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
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

static void show_usage(string program_name, string version_number){
    cout << "Name: Compact Reference Panel Generator V" << version_number << "\n";
    cout << "Usage: " << program_name << " <Option(s)>\n";
    cout << "Option(s):\n";
    cout << "\t-h,--help\t\t\t\t\t\tShow this help message\n";
    cout << "\t-p,--panel <INPUT PANEL FILE>\t\t\t\tInput panel path and file name\n";
    cout << "\t-a,--ancestry <INPUT POPULATION LABEL FILE>\t\tInput panel population label path and file name\n";
    cout << "\t-o,--output <OUTPUT COMPACT PANEL FILE>\t\t\tOutput compact panel path and file name\n";
    cout << "\t-c,--compact <OUTPUT COMPACT POPULATION LABEL FILE>\tOutput compact panel population label path and file name\n";
}

int main(int argc, char *argv[]){
    try {
        cout << "start program" << endl;

        //input variables with default values
        string chromosome_id = "."; //chromosome id
        bool has_chromosome_id = false;

        //other variables
        string version_number = "0.1"; //program version number
        int number_of_populations = 0; //number of ways of admixture population
        char population_delimiter = '\t'; //text format
        char vcf_delimiter = '\t'; //VCF format

        unordered_map<string, int> sample_populations; //key: sample_id, value: population_id
        vector<unordered_map<int, double>> population_zero_frequency_per_site; //population allele frequencies in reference panel; each site has a map having key: population id; value: allele frequency of site having zero value

        int number_of_panel_individuals = 0;
        int number_of_panel_sites = 0;
        vector<string> panel_individual_ids;
        vector<int> panel_physical_positions;
        vector<string> panel_rsids;
        vector<string> panel_ref_alleles;
        vector<string> panel_alt_alleles;

        string input_reference_panel_path_and_file_name; //input reference population panel vcf file
        string input_population_path_and_file_name; //input reference population panel label file
        string output_compact_panel_path_and_file_name; //output reference population panel vcf file
        string output_compact_panel_population_label_path_and_file_name; //output reference population panel label file
        int input_buffer_size = 1;

        //get command line arguments
        string program_name = argv[0];
        string program_arguments = "";
        for (int i = 1; i < argc; i++){
            string argument = argv[i];
            if ((argument == "-h") || (argument == "--help")){
                show_usage(argv[0], version_number);
                cout << "end program" << endl;
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
                    cout << "end program" << endl;
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
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-o") || (argument == "--output")){
                if (i + 1 < argc){
                    output_compact_panel_path_and_file_name = argv[i + 1];
                    program_arguments += " o=" + output_compact_panel_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-o (or --output) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-c") || (argument == "--compact")){
                if (i + 1 < argc){
                    output_compact_panel_population_label_path_and_file_name = argv[i + 1];
                    program_arguments += " c=" + output_compact_panel_population_label_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-c (or --compact) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-b") || (argument == "--buffer")){
                if (i + 1 < argc){
                    input_buffer_size = atoi(argv[i + 1]);
                    program_arguments += " b=" + to_string(input_buffer_size);
                    i++;
                }
                else {
                    cout << "-b (or --buffer) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else {
                cout << "unrecognized arguments: " << argument << endl;
                show_usage(argv[0], version_number);
                cout << "end program" << endl;
                return 1;
            }
        }
        program_arguments = program_arguments.size() > 0 ? program_arguments.substr(1) : "";

        //output variables
        cout << "program name: " << program_name << "\nprogram arguments: " << program_arguments << endl;
        cout << "output files:" << endl;
        cout << "compact reference panel: " << output_compact_panel_path_and_file_name << endl;
        cout << "compact reference panel population label: " << output_compact_panel_population_label_path_and_file_name << endl;
        
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

        cout << "start read reference ancestry data" << endl;

        //use two maps as bidirectional map to store origin labels and ids
        unordered_map<string, int> origin_to_ids;
        unordered_map<int, string> origin_to_labels;
        int label_id = 0;
        set<string> population_labels;

        //read sample population label from file
        input_file_data.open(input_population_path_and_file_name);
        if (!input_file_data){
            cout << "cannot open file " + input_population_path_and_file_name << endl;
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
                        population_labels.insert(population_label);
                    }
                    line_number++;
                }
            }
            input_file_data.close();
        }
        number_of_populations = population_labels.size();
        cout << "number_of_populations=" << number_of_populations << endl;

        cout << "end read reference ancestry data" << endl;

        cout << "start read panel data" << endl;

        //read panel data from file
        streamsize io_buffer_size = input_buffer_size * 1024 * 1024; // default is 1 MB buffer size
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
            cout << "cannot open file " + input_reference_panel_path_and_file_name << endl;
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
                        //get reference SNP cluster id
                        if (i == id_location_on_file){
                            panel_rsids.push_back(tokens[i]);
                        }

                        //get reference allele
                        if (i == reference_base_location_on_file){
                            panel_ref_alleles.push_back(tokens[i]);
                        }

                        //get alternate allele
                        if (i == alternate_base_location_on_file){
                            panel_alt_alleles.push_back(tokens[i]);
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
                        cout << "data at line " + to_string(line_number) + " is not valid" << endl;
                    }
                }
                line_number++;
            }
            input_panel_file.close();
        }
        
        cout << "end read panel data" << endl;

        cout << "number_of_individuals=" << number_of_panel_individuals << endl;
        cout << "number_of_sites=" << number_of_panel_sites << endl;

        //output compact reference panel
        if (chromosome_id.substr(0, 3).compare("chr") || chromosome_id.substr(0, 3).compare("Chr")){
            chromosome_id = chromosome_id.substr(3, chromosome_id.size() - 3);
        }
        output_file_data.open(output_compact_panel_path_and_file_name, ios::trunc);
        if (!output_file_data){
            cout << "cannot create or open file " + output_compact_panel_path_and_file_name << endl;
            exit(1);
        }
        if (output_file_data.is_open()){
            output_file_data << "##fileformat=VCFv4.2" << endl;
            output_file_data << "##source=CompactReferencePopulationPanelv" << version_number << endl;
            output_file_data << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
            output_file_data << "##contig=<ID=" << chromosome_id << ">" << endl;
            output_file_data << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
            for (int i = 0; i < origin_to_labels.size(); i++){
                output_file_data << origin_to_labels[i];
                if (i < origin_to_labels.size() - 1){
                    output_file_data << "\t";
                }
                else {
                    output_file_data << endl;
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
                if (output_file_data.is_open()){
                    output_file_data << chromosome_id << "\t" << panel_physical_positions[i] << "\t" << panel_rsids[i] << "\t" << panel_ref_alleles[i] << "\t" << panel_alt_alleles[i] << "\t.\tPASS\t.\tGT\t";
                    for (int j = 0; j < individuals_site_value.size(); j++){
                        output_file_data << individuals_site_value[j];
                        if (j < individuals_site_value.size() - 1){
                            output_file_data << "\t";
                        }
                        else {
                            output_file_data << endl;
                        }
                    }
                }
            }
        }
        output_file_data.close();

        output_file_data.open(output_compact_panel_population_label_path_and_file_name, ios::trunc);
        if (!output_file_data){
            cout << "cannot create or open file " + output_compact_panel_population_label_path_and_file_name << endl;
            exit(1);
        }
        if (output_file_data.is_open()){
            output_file_data << "#Sample_ID\tAncestry_Label" << endl;
            for (int i = 0; i < origin_to_labels.size(); i++){
                output_file_data << origin_to_labels[i] << "\t" << origin_to_labels[i] << endl;
            }
        }
        output_file_data.close();
        
        cout << "end program" << endl;
    }
    catch (exception & exception_output){
        cerr << exception_output.what() << endl;
    }
    return 0;
}