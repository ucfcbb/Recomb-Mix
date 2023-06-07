//  * --------------------------------------------------------------------------------------------------------
//  * Name: RecombMix.cpp
//  * Description: Local Ancestry Inference based on improved Loter with new model, using recombination rate.
//  * Author: Yuan Wei 
//  * Created on: Jan 21, 2023
//  * Modified on: Jun 02, 2023
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

static void show_usage(string program_name){
    cout << "Usage: " << program_name << " <Option(s)>\n";
    cout << "Option(s):\n";
    cout << "\t-h,--help\t\t\t\t\tShow this help message\n";
    cout << "\t-p,--panel <INPUT PANEL FILE>\t\t\tInput panel path and file name\n";
    cout << "\t-q,--query <INPUT QUERY FILE>\t\t\tInput query path and file name\n";
    cout << "\t-g,--genetic <INPUT GENETIC MAPPING FILE>\tInput genetic mapping path and file name\n";
    cout << "\t-a,--ancestry <INPUT POPULATION ANCESTRY FILE>\tInput population ancestry path and file name\n";
    cout << "\t-o,--output <OUTPUT DIRECTORY PATH>\t\tOutput directory path for all files\n";
    cout << "\t-e,--weight <WEIGHT>\t\t\t\tWeight of recombination rate in cost function\n";
    cout << "\t-u,--outputcompactpanel <IDENTIFIER>\t\tSpecify whether output compact reference panel (1: true; 0: false)\n";
}

int main(int argc, char *argv[]){
    try {
        cout << "start program" << endl;

        //input variables with default values
        string chromosome_id = "."; //chromosome id
        bool has_chromosome_id = false;
        string output_directory_path = "./"; //default is current folder
        double weight = 1.5; //recombination rate weight in cost function

        //other variables
        string version_number = "0.4"; //program version number
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

        vector<vector<bool>> query_sites;
        int number_of_query_individuals = 0;
        int number_of_query_sites = 0;
        vector<string> query_individual_ids;
        vector<int> query_physical_positions;

        vector<vector<bool>> panel_sites;
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
        bool output_compact_panel = false;
        bool has_queries = false;

        //get command line arguments
        string program_name = argv[0];
        string program_arguments = "";
        for (int i = 1; i < argc; i++){
            string argument = argv[i];
            if ((argument == "-h") || (argument == "--help")){
                show_usage(argv[0]);
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
            else if ((argument == "-q") || (argument == "--query")){
                if (i + 1 < argc){
                    input_query_panel_path_and_file_name = argv[i + 1];
                    has_queries = true;
                    program_arguments += " q=" + input_query_panel_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-q (or --query) option requires one argument\n";
                    cout << "end program" << endl;
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
                    output_directory_path = argv[i + 1];
                    program_arguments += " o=" + output_directory_path;
                    i++;
                }
                else {
                    cout << "-o (or --output) option requires one argument\n";
                    cout << "end program" << endl;
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
                    cout << "end program" << endl;
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
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else {
                cout << "unrecognized arguments: " << argument << endl;
                show_usage(argv[0]);
                cout << "end program" << endl;
                return 1;
            }
        }
        program_arguments = program_arguments.size() > 0 ? program_arguments.substr(1) : "";

        //construct output files
        if (output_directory_path.substr(output_directory_path.size() - 1).compare("/") == 0){
            output_directory_path = output_directory_path.substr(0, output_directory_path.size() - 1);
        }
        output_inference_individuals_path_and_file_name = output_directory_path + "/admix_inferred_ancestral_values_local.txt";
        output_compact_panel_path_and_file_name = output_directory_path + "/compact_reference_panel.vcf.gz";
        output_compact_panel_population_label_path_and_file_name = output_directory_path + "/compact_reference_panel_population_label.txt";

        //output variables
        cout << "program name: " << program_name << "\nprogram arguments: " << program_arguments << endl;
        cout << "weight=" << fixed << setprecision(6) << weight << endl;
        cout << "output files:" << endl;
        cout << "local ancestray inference: " << output_inference_individuals_path_and_file_name << endl;
        if (output_compact_panel){
            cout << "compact reference panel: " << output_compact_panel_path_and_file_name << endl;
            cout << "compact reference panel population label: " << output_compact_panel_population_label_path_and_file_name << endl;
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
                        sample_index_populations[line_number] = origin_to_ids[population_label];
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
        stream_in_buffer_panel.push(input_panel_file);
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
                        //store current site values from panel individuals
                        panel_sites.push_back(site_of_panel_individuals);
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

        //output compact reference panel
        if (output_compact_panel){
            ofstream output_panel_file(output_compact_panel_path_and_file_name, ios_base::out | ios_base::binary);
            if (!output_panel_file){
                cout << "cannot create or open file " + output_compact_panel_path_and_file_name << endl;
                exit(1);
            }
            boost::iostreams::filtering_streambuf<boost::iostreams::output> stream_out_buffer;
            stream_out_buffer.push(boost::iostreams::gzip_compressor());
            stream_out_buffer.push(output_panel_file);
            ostream output_vcf_file_data(&stream_out_buffer);
            if (output_panel_file.is_open()){
                output_vcf_file_data << "##fileformat=VCFv4.2" << endl;
                output_vcf_file_data << "##source=RecombMix_v" << version_number << endl;
                output_vcf_file_data << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
                output_vcf_file_data << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
                for (int i = 0; i < origin_to_labels.size(); i++){
                    output_vcf_file_data << origin_to_labels[i];
                    if (i < origin_to_labels.size() - 1){
                        output_vcf_file_data << "\t";
                    }
                    else {
                        output_vcf_file_data << endl;
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
                    output_vcf_file_data << chromosome_id << "\t" << panel_physical_positions[i] << "\t.\tA\tC\t.\tPASS\t.\tGT\t";
                    for (int j = 0; j < individuals_site_value.size(); j++){
                        output_vcf_file_data << individuals_site_value[j];
                        if (j < individuals_site_value.size() - 1){
                            output_vcf_file_data << "\t";
                        }
                        else {
                            output_vcf_file_data << endl;
                        }
                    }
                }
            }
            boost::iostreams::close(stream_out_buffer);
            output_panel_file.close();

            output_file_data.open(output_compact_panel_population_label_path_and_file_name, ios::trunc);
            if (!output_file_data){
                cout << "cannot create or open file " + output_compact_panel_population_label_path_and_file_name << endl;
                exit(1);
            }
            if (output_file_data.is_open()){
                output_file_data << "#Sample_ID\tPopulation_Label" << endl;
                for (int i = 0; i < origin_to_labels.size(); i++){
                    output_file_data << origin_to_labels[i] << "\t" << origin_to_labels[i] << endl;
                }
            }
            output_file_data.close();
        }

        if (!has_queries){
            cout << "end program" << endl;
            return 0;
        }

        cout << "start read query data" << endl;

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
            cout << "cannot open file " + input_query_panel_path_and_file_name << endl;
            exit(1);
        }
        stream_in_buffer_query.push(input_query_file);
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
                    bool is_a_valid_line = false;
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
                            }
                            else {
                                //this site from query is not in the reference panel, skip it
                                number_of_query_sites++;
                                is_a_valid_line = true;
                                break;
                            }
                        }
                        //get query individual first and second haplotype site values
                        if (tokens[i].substr(1, 1).compare("|") == 0){ //phased
                            site_of_query_individuals.push_back(stoi(tokens[i].substr(0, 1)));
                            site_of_query_individuals.push_back(stoi(tokens[i].substr(2, 1)));
                        }
                    }
                    if (site_of_query_individuals.size() > 0){
                        //store current site values from query individuals
                        query_sites.push_back(site_of_query_individuals);
                        site_of_query_individuals.clear();
                        number_of_query_sites++;
                    }
                    else {
                        if (!is_a_valid_line){
                            cout << "data at line " + to_string(line_number) + " is not valid" << endl;
                        }
                    }
                }
                line_number++;
            }
            input_query_file.close();
        }

        cout << "end read query data" << endl;

        cout << "number_of_intersecting_markers=" << panel_physical_position_index_filtered.size();
        cout << ", percentage_of_markers_in_query_covered_by_reference=" << fixed << setprecision(2) << (double(panel_physical_position_index_filtered.size()) / double(query_physical_positions.size())) * 100 << "%" << endl;

        cout << "start read map data" << endl;

        //read genetic map data from file
        input_file_data.open(input_map_path_and_file_name);
        if (!input_file_data){
            cout << "cannot open file " + input_map_path_and_file_name << endl;
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
        double rr_min = numeric_limits<double>::max();
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

        cout << "end read map data" << endl;

        cout << "start verify data" << endl;

        bool is_data_satisfied = true;

        //verify parameters
        if (panel_physical_position_index_filtered.size() <= 0){
            is_data_satisfied = false;
            cout << "reference panel has no intersecting markers with query" << endl;
        }
        if (number_of_populations == 0){
            is_data_satisfied = false;
            cout << "reference panel has no ancestry information" << endl;
        }
        if (weight < 0.0){
            is_data_satisfied = false;
            cout << "weight of recombination rate is invalid" << endl;
        }

        cout << "end verify data" << endl;

        if (!is_data_satisfied){
            cout << "end program" << endl;
            return 0;
        }

        int number_of_query_haplotypes = number_of_query_individuals * 2;
        int number_of_panel_haplotypes = number_of_panel_individuals * 2;
        cout << "number of query haplotypes: " << to_string(number_of_query_haplotypes) << endl;
        cout << "number of query sites: " << to_string(number_of_query_sites) << endl;
        cout << "number of panel haplotypes: " << to_string(number_of_panel_haplotypes) << endl;
        cout << "number of panel sites: " << to_string(number_of_panel_sites) << endl;

        int m = number_of_panel_haplotypes;
        int n = panel_physical_position_index_filtered.size();
        int q = number_of_query_haplotypes;

        //initialize graph
        vector<vector<double>> scores; //stores score of each site
        vector<unordered_map<int, tuple<double, int>>> scores_min; //stores minimum scores and paths of each population appearing on the site
        vector<vector<int>> paths; //store optimal path of each site (haplotype index of the previous site)

        double score_final = numeric_limits<double>::max();
        int last_site_haplotype_index = -1;
        vector<int> path_backtracks;

        cout << "start run queries" << endl;

        output_file_data.open(output_inference_individuals_path_and_file_name, ios::app);
        if (!output_file_data){
            cout << "cannot create or open file " + output_inference_individuals_path_and_file_name << endl;
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
                    output_file_data << endl;
                }
            }
        }
        output_file_data.close();
        
        for (int query_haplotype_index = 0; query_haplotype_index < q; query_haplotype_index++){
            //get query haplotype
            vector<bool> query_site;
            for (int j = 0 ; j < n; j++){
                query_site.push_back(query_sites[j][query_haplotype_index]);
            }

            //build compact graph based on compact reference panel
            for (int i = 0; i < panel_physical_position_index_filtered.size(); i++){
                int node_site_index = i;
                int prev_node_site_index = i > 0 ? i - 1 : 0;
                double recombination_rate_penalty = panel_recombination_rate_penalties[prev_node_site_index];
                vector<double> hss(m, numeric_limits<double>::max()); //haplotype scores of the site
                unordered_map<int, tuple<double, int>> hss_min; //min haplotype scores paths of the site per population; key: population_label_id, value: score, haplotype_index
                vector<int> hps(m, -1); //optimal haplotype paths of the site
                int population_label_id = -1; //"UNKNOWN";
                for (const auto & [key_population_label_id, value_allele_frequency_zero] : population_zero_frequency_per_site[panel_physical_position_index_filtered[i]]){
                    population_label_id = key_population_label_id;
                    double score = numeric_limits<double>::max();
                    int path_index = -1;
                    double score_l = 0.0;
                    int haplotype_id_of_site = -1;
                    //check site having allele value zero and one
                    if (value_allele_frequency_zero > 0.0){
                        int panel_site_value = 0;
                        int haplotype_id_of_site = population_label_id * 2 + 0;
                        double score = numeric_limits<double>::max();
                        int path_index = -1;
                        double score_l = 0.0;
                        if (query_site[node_site_index] == panel_site_value){
                            score_l = 0.0;
                        }
                        else {
                            score_l = 1.0;
                        }
                        if (i == 0){
                            score = score_l;
                        }
                        else {
                            for (auto const & [population_label_id_prev, score_path] : scores_min[i - 1]){
                                double score_r = population_label_id_prev == population_label_id ? 0.0 : recombination_rate_penalty;
                                score_r = (weight * score_r);
                                if (score > get<0>(score_path) + score_l + score_r){
                                    score = get<0>(score_path) + score_l + score_r;
                                    path_index = get<1>(score_path);
                                }
                            }
                        }

                        //update trackers
                        hss[haplotype_id_of_site] = score;
                        if (hss_min.size() > 0 && hss_min.count(population_label_id) > 0){
                            if (get<0>(hss_min[population_label_id]) > score){
                                hss_min[population_label_id] = tuple<double, int>(score, haplotype_id_of_site);
                            }
                        }
                        else {
                            hss_min[population_label_id] = tuple<double, int>(score, haplotype_id_of_site);
                        }
                        hps[haplotype_id_of_site] = path_index;

                        //track last score and path
                        if (i == panel_physical_position_index_filtered.size() - 1){
                            if (score_final > score){
                                score_final = score;
                                last_site_haplotype_index = haplotype_id_of_site;
                            }
                        }
                    }
                    if (1.0 - value_allele_frequency_zero > 0.0){
                        int panel_site_value = 1;
                        int haplotype_id_of_site = population_label_id * 2 + 1;
                        double score = numeric_limits<double>::max();
                        int path_index = -1;
                        double score_l = 0.0;
                        if (query_site[node_site_index] == panel_site_value){
                            score_l = 0.0;
                        }
                        else {
                            score_l = 1.0;
                        }
                        if (i == 0){
                            score = score_l;
                        }
                        else {
                            for (auto const & [population_label_id_prev, score_path] : scores_min[i - 1]){
                                double score_r = population_label_id_prev == population_label_id ? 0.0 : recombination_rate_penalty;
                                score_r = (weight * score_r);
                                if (score > get<0>(score_path) + score_l + score_r){
                                    score = get<0>(score_path) + score_l + score_r;
                                    path_index = get<1>(score_path);
                                }
                            }
                        }

                        //update trackers
                        hss[haplotype_id_of_site] = score;
                        if (hss_min.size() > 0 && hss_min.count(population_label_id) > 0){
                            if (get<0>(hss_min[population_label_id]) > score){
                                hss_min[population_label_id] = tuple<double, int>(score, haplotype_id_of_site);
                            }
                        }
                        else {
                            hss_min[population_label_id] = tuple<double, int>(score, haplotype_id_of_site);
                        }
                        hps[haplotype_id_of_site] = path_index;

                        //track last score and path
                        if (i == panel_physical_position_index_filtered.size() - 1){
                            if (score_final > score){
                                score_final = score;
                                last_site_haplotype_index = haplotype_id_of_site;
                            }
                        }
                    }
                }
                scores.push_back(hss);
                scores_min.push_back(hss_min);
                paths.push_back(hps);
            }

            int j = last_site_haplotype_index;
            for (int i = scores.size() - 1; i >= 0; i--){
                path_backtracks.push_back(j);
                j = paths[i][j];
            }

            //merge and store consecutive sites having the same ancestral value
            vector<tuple<int, int, int>> inferred_segments; //inferred ancestral value of each site (compressed): start_index, end_index (inclusive), population_label_id
            int start_position_prev = -2;
            int population_label_prev = -2;
            int k1 = query_physical_position_index_filtered.size() - 1; //for tracking intersected sites
            int k2 = panel_physical_position_index_filtered.size() - 1; //for getting population label
            for (int i = number_of_query_sites - 1; i >= 0; i--){ //traverse original number of sites in query panel; if no intersection with reference panel, mark -1 as its population label
                //get population label
                int population_label_id_selected = -1;
                if (query_physical_position_index_filtered[k1] == i){
                    //intersected
                    while (panel_physical_position_index_filtered[k2] > i){
                        //it is guaranteed to find a site which satisfied panel_physical_position_index_filtered[k2] == i
                        k2--;
                    }
                    population_label_id_selected = path_backtracks[k2] / 2;
                    k1--;
                }
                else {
                    //not intersected
                    population_label_id_selected = -1;
                }
                //merge consecutive sites having the same ancestral value
                if (population_label_id_selected != population_label_prev){
                    if (population_label_prev != -2){
                        //store the segment inferred result
                        inferred_segments.push_back(tuple<int, int, int>(start_position_prev, query_physical_positions[number_of_query_sites - 1 - i - 1], population_label_prev));
                    }
                    start_position_prev = query_physical_positions[number_of_query_sites - 1 - i];
                    population_label_prev = population_label_id_selected;
                }
            }
            //store the last segment inferred result
            inferred_segments.push_back(tuple<int, int, int>(start_position_prev, query_physical_positions[number_of_query_sites - 1], population_label_prev));

            //output inference result
            string query_individual_id = query_individual_ids[query_haplotype_index / 2];
            output_file_data.open(output_inference_individuals_path_and_file_name, ios::app);
            if (!output_file_data){
                cout << "cannot create or open file " + output_inference_individuals_path_and_file_name << endl;
                exit(1);
            }
            if (output_file_data.is_open()){
                //output file format: query_individual_id_query_haplotype_index \t origin_label_ids of each segment's start position, end position, population value (merged if consecutive segments having the same ancestral value) separated by \t
                output_file_data << query_individual_id << "_" << (query_haplotype_index % 2);
                for (int i = 0; i < inferred_segments.size(); i++){
                    output_file_data << "\t" << get<0>(inferred_segments[i]) << "\t" << get<1>(inferred_segments[i]) << "\t" << get<2>(inferred_segments[i]);
                }
                output_file_data << endl;
            }
            output_file_data.close();

            //reset variables for next query haplotype
            scores.clear();
            scores_min.clear();
            paths.clear();
            score_final = numeric_limits<double>::max();
            last_site_haplotype_index = -1;
            path_backtracks.clear();
        }

        cout << "end run queries" << endl;

        cout << "end program" << endl;
    }
    catch (exception & exception_output){
        cerr << exception_output.what() << endl;
    }
    return 0;
}