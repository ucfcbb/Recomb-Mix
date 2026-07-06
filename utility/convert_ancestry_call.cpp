//  * --------------------------------------------------------------------------------------------------------
//  * Name: convert_ancestry_call.cpp
//  * Description: Convert local ancestry inference result file format from run-length encoding to ancestry state or ancestry dosage.
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

static void show_usage(string program_name, string version_number){
    cout << "Name: Ancestry Call Converter V" << version_number << "\n";
    cout << "Usage: " << program_name << " <Option(s)>\n";
    cout << "Option(s):\n";
    cout << "\t-h,--help\t\t\t\t\tShow this help message\n";
    cout << "\t-i,--input <INPUT RESULT FILE>\t\t\tInput local ancestry inference result path and file name (run-length encoding format)\n";
    cout << "\t-o,--output <OUTPUT RESULT FILE>\t\tOutput local ancestry inference result path and file name\n";
    cout << "\t-f,--format <IDENTIFIER>\t\t\tOutput local ancestry inference result file format (0: ancestry state; 1: ancestry dosage)\n";
}

int main(int argc, char *argv[]){
    try {
        cout << "start program" << "\n";

        //input variables with default values
        string version_number = "0.1"; //program version number
        int number_of_populations = 0; //number of ways of admixture population
        int number_of_query_haplotypes = 0;
        char file_delimiter = '\t'; //local ancestry inference result file format

        int number_of_query_individuals = 0;
        int number_of_query_sites = 0;
        vector<string> query_individual_ids; //individual haplotype ids
        vector<string> query_individual_names; //individual ids
        vector<int> query_physical_positions;

        char inference_delimiter = '\t';
        string input_inference_path_and_file_name;
        string output_inference_path_and_file_name = "./inferred_local_ancestral_values.tsv.gz"; //default path and file name

        int input_buffer_size = 1;
        int output_format = 0;

        //get command line arguments
        string program_name = argv[0];
        string program_arguments = "";
        for (int i = 1; i < argc; i++){
            string argument = argv[i];
            if ((argument == "-h") || (argument == "--help")){
                show_usage(argv[0], version_number);
                cout << "end program" << "\n";
                return 0;
            }
            else if ((argument == "-i") || (argument == "--input")){
                if (i + 1 < argc){
                    input_inference_path_and_file_name = argv[i + 1];
                    program_arguments += " i=" + input_inference_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-i (or --input) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-o") || (argument == "--output")){
                if (i + 1 < argc){
                    output_inference_path_and_file_name = argv[i + 1];
                    program_arguments += " o=" + output_inference_path_and_file_name;
                    i++;
                }
                else {
                    cout << "-o (or --output) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else if ((argument == "-f") || (argument == "--format")){
                if (i + 1 < argc){
                    output_format = stoi(argv[i + 1]);
                    program_arguments += " f=" + to_string(output_format);
                    i++;
                }
                else {
                    cout << "-f (or --format) option requires one argument\n";
                    cout << "end program" << "\n";
                    return 1;
                }
            }
            else {
                cout << "unrecognized arguments: " << argument << "\n";
                show_usage(argv[0], version_number);
                cout << "end program" << "\n";
                return 1;
            }
        }
        program_arguments = program_arguments.size() > 0 ? program_arguments.substr(1) : "";

        //output variables
        cout << "program name: " << program_name << "\nprogram arguments: " << program_arguments << "\n";
        cout << "output file format: ";
        if (output_format == 1){
            cout << "ancestry dosage" << endl;
        }
        else if (output_format == 0){
            cout << "ancestry state" << endl;
        }
        else {
            cout << "not supported" << endl;
        }
        if (output_inference_path_and_file_name.substr(output_inference_path_and_file_name.size() - 7, 7).compare(".tsv.gz") != 0){
            output_inference_path_and_file_name += ".tsv.gz";
        }
        cout << "output file name: ";
        cout << output_inference_path_and_file_name << endl;
        
        //variables for input and output
        ifstream input_file_data;
        ofstream output_file_data;

        cout << "start read inference data" << "\n";
        //read inferred data from file (run-length encoding format)
        map<int, string> ancestry_labels; //key: ancestry id, value: ancestry label
        vector<int> physical_positions; //panel physical positions
        unordered_map<int, int> physical_positions_site_index; //key: physical_position, value: site_index
        vector<string> query_output_buffer;
        int section_id = 0;
        input_file_data.open(input_inference_path_and_file_name);
        if (!input_file_data){
            cout << "cannot open file " + input_inference_path_and_file_name << endl;
            exit(1);
        }
        if (input_file_data.is_open()){
            int line_number = 0;
            string line_from_file;
            while (getline(input_file_data, line_from_file)){
                if (line_from_file.substr(0, 9).compare("#Ancestry") == 0){
                    section_id = 0;
                }
                else if (line_from_file.substr(0, 9).compare("#Position") == 0){
                    section_id = 1;
                }
                else if (line_from_file.substr(0, 7).compare("#Result") == 0){
                    section_id = 2;
                }
                else {
                    //read data
                    if (section_id == 0){
                        //read ancestry section
                        vector<string> tokens = split_string(line_from_file, inference_delimiter);
                        if (tokens.size() >= 2){
                            //format: ancestry label \t ancestry id
                            string ancestry_label = remove_spaces(tokens[0]);
                            int ancestry_id = stoi(tokens[1]);
                            ancestry_labels[ancestry_id] = ancestry_label;
                        }
                    }
                    else if (section_id == 1){
                        //read position section
                        int physical_location = stoi(line_from_file);
                        int physical_position_site_index = physical_positions.size();
                        physical_positions.push_back(physical_location);
                        physical_positions_site_index[physical_location] = physical_position_site_index;
                    }
                    else {
                        //read result section
                        query_output_buffer.push_back(line_from_file);
                    }
                    line_number++;
                }
            }
            input_file_data.close();
        }

        cout << "end read inference data" << endl;

        cout << "start output inference data" << endl;

        number_of_populations = ancestry_labels.size();
        number_of_query_haplotypes = query_output_buffer.size();

        if (output_format == 0){
            //convert output format from run-length encoding to ancestry state
            vector<vector<int>> ancestry_state(physical_positions_site_index.size(), vector<int>(number_of_populations * number_of_query_haplotypes, 0));
            for (int query_haplotype_index = 0; query_haplotype_index < number_of_query_haplotypes; query_haplotype_index++){   
                vector<string> tokens = split_string(query_output_buffer[query_haplotype_index], inference_delimiter);
                query_individual_ids.push_back(tokens[0]);
                int token_index = 1;
                while (token_index + 2 < tokens.size()){
                    int start_physical_position = stoi(tokens[token_index]);
                    int end_physical_position = stoi(tokens[token_index + 1]);
                    int ancestry_value = stoi(tokens[token_index + 2]);
                    int start_physical_position_index = physical_positions_site_index[start_physical_position];
                    int end_physical_position_index = physical_positions_site_index[end_physical_position];
                    for (int i = start_physical_position_index; i <= end_physical_position_index; i++){
                        ancestry_state[i][number_of_populations * query_haplotype_index + ancestry_value] = 1;
                    }
                    token_index += 3;
                }
            }

            //output in ancestry state format
            streamsize io_buffer_size = input_buffer_size * 1024 * 1024; // default is 1 MB buffer size
            ofstream output_result_file(output_inference_path_and_file_name, ios_base::out | ios_base::binary);
            if (!output_result_file){
                cout << "cannot create or open file " + output_inference_path_and_file_name << "\n";
                exit(1);
            }
            boost::iostreams::filtering_streambuf<boost::iostreams::output> stream_out_buffer;
            stream_out_buffer.push(boost::iostreams::gzip_compressor());
            stream_out_buffer.push(output_result_file, io_buffer_size);
            ostream output_result_file_data(&stream_out_buffer);
            if (output_result_file.is_open()){
                output_result_file_data << "Position";
                for (int i = 0; i < query_individual_ids.size(); i++){
                    for (const auto & [id, lbl] : ancestry_labels){
                        output_result_file_data << "\t" << query_individual_ids[i] + "_" + lbl;
                    }
                }
                output_result_file_data << "\n";
                for (int i = 0; i < ancestry_state.size(); i++){
                    output_result_file_data << physical_positions[i];
                    for (int j = 0; j < ancestry_state[i].size(); j++){
                        output_result_file_data << "\t" << ancestry_state[i][j];
                    }
                    output_result_file_data << "\n";
                }
            }
            boost::iostreams::close(stream_out_buffer);
            output_result_file.close();
        }
        else if (output_format == 1){
            //convert output format from run-length encoding to ancestry dosage
            vector<vector<int>> anc_dos(physical_positions_site_index.size(), vector<int>(number_of_populations * (number_of_query_haplotypes / 2), 0));
            for (int query_haplotype_index = 0; query_haplotype_index < number_of_query_haplotypes; query_haplotype_index++){   
                vector<string> tokens = split_string(query_output_buffer[query_haplotype_index], inference_delimiter);
                string individual_id_with_haplotype_id = tokens[0];
                size_t last_position = individual_id_with_haplotype_id.rfind('_');
                string individual_name = individual_id_with_haplotype_id.substr(0, last_position); //excluding the '_' character
                if (query_individual_names.size() == 0 || query_individual_names[query_individual_names.size() - 1] != individual_name){
                    query_individual_names.push_back(individual_name);
                }
                int token_index = 1;
                while (token_index + 2 < tokens.size()){
                    int start_physical_position = stoi(tokens[token_index]);
                    int end_physical_position = stoi(tokens[token_index + 1]);
                    int ancestry_value = stoi(tokens[token_index + 2]);
                    int start_physical_position_index = physical_positions_site_index[start_physical_position];
                    int end_physical_position_index = physical_positions_site_index[end_physical_position];
                    for (int i = start_physical_position_index; i <= end_physical_position_index; i++){
                        anc_dos[i][number_of_populations * (query_haplotype_index / 2) + ancestry_value] += 1;
                    }
                    token_index += 3;
                }
            }
            //output in ancestry dosage format
            streamsize io_buffer_size = input_buffer_size * 1024 * 1024; // default is 1 MB buffer size
            ofstream output_result_file(output_inference_path_and_file_name, ios_base::out | ios_base::binary);
            if (!output_result_file){
                cout << "cannot create or open file " + output_inference_path_and_file_name << "\n";
                exit(1);
            }
            boost::iostreams::filtering_streambuf<boost::iostreams::output> stream_out_buffer;
            stream_out_buffer.push(boost::iostreams::gzip_compressor());
            stream_out_buffer.push(output_result_file, io_buffer_size);
            ostream output_result_file_data(&stream_out_buffer);

            if (output_result_file.is_open()){
                output_result_file_data << "Position";
                for (int i = 0; i < query_individual_names.size(); i++){
                    for (const auto & [id, lbl] : ancestry_labels){
                        output_result_file_data << "\t" << query_individual_names[i] + "_" + lbl;
                    }
                }
                output_result_file_data << "\n";
                for (int i = 0; i < anc_dos.size(); i++){
                    output_result_file_data << physical_positions[i];
                    for (int j = 0; j < anc_dos[i].size(); j++){
                        output_result_file_data << "\t" << anc_dos[i][j];
                    }
                    output_result_file_data << "\n";
                }
            }
            boost::iostreams::close(stream_out_buffer);
            output_result_file.close();
        }
        else {
            cout << "format " << output_format << " not supported" << endl;
        }

        cout << "end output inference data" << endl;

        cout << "end program" << "\n";
    }
    catch (exception & exception_output){
        cerr << exception_output.what() << "\n";
    }
    return 0;
}