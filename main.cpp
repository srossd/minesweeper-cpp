#include <stdio.h>
#include <random>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <queue>
#include <bitset>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <memory>
#include <stdexcept>

using namespace std;

typedef struct {
    pair<int, int> coord;
    bool is_mine;
    int core;
    bool success;
} inference;

vector<vector<int>> buildGrid(int, double, bool = false);
void pretty_print(vector<vector<int>>&);
vector<vector<int>> eff_grid(vector<vector<int>>&);
void flag(vector<vector<int>>&, int, int);
vector<pair<int,int>> frontier_inner(vector<vector<int>>&);
vector<pair<int,int>> frontier_outer(vector<vector<int>>&);
inference infer(vector<vector<int>>&);
pair<bool, int> inference_play(vector<vector<int>>&, bool=false, bool=false);
vector<pair<int,int>> neighbors(int, int, int);
void reveal(vector<vector<int>>&, vector<vector<int>>&, int, int);
void test_muser();
void get_data(int, int);
string exec(const char*);
pair<int, bool> call_muser(vector<pair<int, vector<int>>>, int);

int main() {
    get_data(20, 100);
    get_data(40, 100);
    get_data(10, 100);
}

pair<int, bool> call_muser(vector<pair<int, vector<int>>> clauses, int nvars) {
    ofstream gcnf("temp.gcnf");
    gcnf << "p gcnf " << nvars << " " << clauses.size() << " " << clauses.back().first << endl;
    for(pair<int,vector<int>> clause : clauses) {
        gcnf << "{" << clause.first << "} ";
        for(int var : clause.second)
            gcnf << var << " ";
        gcnf << "0" << endl;
    }
    gcnf.close();

    string muser_output = exec("./muser2 -grp -test temp.gcnf");
    size_t result_idx = muser_output.find("result:");
    bool sat = muser_output[result_idx+8] == 'S';

    size_t size_idx = muser_output.find("MUS size:");
    int size;
    sscanf(muser_output.substr(size_idx + 10).c_str(), "%d", &size);

    return pair<int,bool>(size, sat);
}

string exec(const char* cmd) {
    array<char, 128> buffer;
    string result;
    unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

void get_data(int N, int trials) {
    vector<double> ps(30);
    for(int i = 0; i < ps.size(); i++)
        ps[i] = 0.01 + (0.5-0.01)*i/(ps.size()-1);

    
    vector<double> data(ps.size(), 0);
    vector<double> timing(ps.size(), 0);
    vector<double> avg_max_core(ps.size(), 0);

    printf("N = %d\n", N);
    for(int j = 0; j < ps.size(); j++) {
        double p = ps[j];
        printf("p = %.2f (%d/%d)\n",p, (j+1), ps.size());

        time_t start = time(NULL);
        for(int i = 0; i < trials; i++) {
            vector<vector<int>> grid = buildGrid(N, p, true);
            pair<bool,int> won = inference_play(grid);
            if(won.first)
                data[j]++;
            avg_max_core[j] += won.second;
        }
        timing[j] = (time(NULL)-start+0.0)/trials;
        data[j] /= trials;
        avg_max_core[j] /= trials;
        printf("Solved %0.0f%% of games, took %0.0f seconds, average max core = %0.2f\n\n", data[j]*100, timing[j]*trials, avg_max_core[j]);
    }

    ofstream data_file("data/solvability_data_"+to_string(N)+".txt");
    for (const auto &e : data) data_file << e << "\n";

    ofstream timing_file("data/timing_data_"+to_string(N)+".txt");
    for (const auto &e : timing) timing_file << e << "\n";

    ofstream core_file("data/core_data_"+to_string(N)+".txt");
    for (const auto &e : avg_max_core) core_file << e << "\n";
}

int count_mines(vector<vector<int>> grid) {
    int ans = 0;
    for(int i = 0; i < grid.size(); i++)
        for(int j = 0; j < grid.size(); j++)
            ans += grid[i][j] == -1;
    return ans;
}

pair<bool, int> inference_play(vector<vector<int>>& grid, bool display, bool verbose) {
    int num_mines = count_mines(grid);
    vector<vector<int>> game(grid.size());
    for(int i = 0; i < grid.size(); i++)
        game[i] = vector<int>(grid[0].size(), -2);
    
    int mc = 0;
    reveal(grid, game, 0, 0);
    while(count_mines(game) < num_mines) {
        if(display)
            pretty_print(game);
        
        inference inf = infer(game);
        if(inf.success) {
            if(inf.is_mine) {
                flag(game, inf.coord.first, inf.coord.second);
                if(verbose)
                    printf("A mine is found at (%d, %d)\n", inf.coord.first, inf.coord.second);
            }
            else {
                reveal(grid, game, inf.coord.first, inf.coord.second);
                if(verbose)
                    printf("A safe square is found at (%d, %d)\n", inf.coord.first, inf.coord.second);
            }
            if(inf.core > mc)
                mc = inf.core;
        }
        else {
            int remaining = num_mines - count_mines(game);
            if(verbose)
                printf("No more inferences are possible. %d mines remain hidden\n", remaining);
            return pair<bool, int>(false, mc);
        }
    }
    if(verbose)
        printf("The game is completed by inference\n");
    return pair<bool, int>(true, mc);
}

inference infer(vector<vector<int>>& game) {
    vector<vector<int>> eg = eff_grid(game);
    vector<pair<int,int>> f_out = frontier_outer(game);
    vector<pair<int,int>> f_in = frontier_inner(game);

    map<pair<int,int>, int> ind_f_out;
    for(int i = 0; i < f_out.size(); i++)
        ind_f_out[f_out[i]] = i+1;

    vector<pair<int,vector<int>>> clauses;
    for(int i = 0; i < f_in.size(); i++) {
        vector<int> bool_vars;
        for(pair<int,int> n : neighbors(f_in[i].first, f_in[i].second, game.size())) {
            if(eg[n.first][n.second] == -2) {
                bool_vars.push_back(ind_f_out[n]);
            }
        }

        if(eg[f_in[i].first][f_in[i].second] == 0) {
            for(pair<int,int> n : neighbors(f_in[i].first, f_in[i].second, game.size())) {
                if(eg[n.first][n.second] == -2) {
                    return {.coord = n, .is_mine = false, .core = 1, .success = true};
                }
            }
        }
        if(eg[f_in[i].first][f_in[i].second] == bool_vars.size()) {
            for(pair<int,int> n : neighbors(f_in[i].first, f_in[i].second, game.size())) {
                if(eg[n.first][n.second] == -2) {
                    return {.coord = n, .is_mine = true, .core = 1, .success = true};
                }
            }
        }

        for(int b = 0; b < 1 << bool_vars.size(); b++) {
            bitset<8> bs(b);
            if(bs.count() != eg[f_in[i].first][f_in[i].second]) {
                vector<int> clause;
                for(int i = 0; i < bool_vars.size(); i++)
                    clause.push_back((1-2*bs[i])*bool_vars[i]);

                clauses.push_back(pair<int,vector<int>>(i, clause));
            }
        } 
    }

    int extra_group = f_in.size();
    vector<int> is_mine_clause;
    pair<int,bool> res;
    for(int i = 0; i < f_out.size(); i++) {
        is_mine_clause = {i+1};
        clauses.push_back(pair<int,vector<int>>(extra_group, is_mine_clause));
        res = call_muser(clauses, f_out.size());
        if(!res.second) {
            return {.coord = f_out[i], .is_mine = false, .core = res.first, .success = true};
        }
        clauses.pop_back();
        
        is_mine_clause = {-(i+1)};
        clauses.push_back(pair<int,vector<int>>(extra_group, is_mine_clause));
        res = call_muser(clauses, f_out.size());
        if(!res.second) {
            return {.coord = f_out[i], .is_mine = true, .core = res.first, .success = true};
        }
        clauses.pop_back();

    }
    
    return {.coord = f_out[0], .is_mine = true, .core = 0, .success = false};
}

void reveal(vector<vector<int>>& grid, vector<vector<int>>& game, int i, int j) {
    queue<pair<int,int>> to_click;
    to_click.push(pair<int,int>(i,j));

    set<pair<int,int>> seen;
    while(!to_click.empty()) {
        pair<int,int> coord = to_click.front();
        to_click.pop();
        if(grid[coord.first][coord.second] != -1) {
            game[coord.first][coord.second] = grid[coord.first][coord.second];

            if(game[coord.first][coord.second] == 0) {
                for(pair<int,int> n : neighbors(coord.first, coord.second, grid.size())) {
                    if(!seen.count(n)) {
                        to_click.push(n);
                        seen.insert(n);
                    }
                }
            }
        }
    }
}

void flag(vector<vector<int>>& game, int i, int j) {
    game[i][j] = -1;
}

vector<pair<int,int>> neighbors(int i, int j, int size) {
    vector<pair<int,int>> n;
    for(int dx = -1; dx <= 1; dx++)
        for(int dy = -1; dy <= 1; dy++)
            if(dx != 0 || dy != 0)
                n.push_back(pair<int,int>((i+dx+size) % size, (j+dy+size) % size));
    return n;
}

vector<pair<int,int>> frontier_inner(vector<vector<int>>& game) {
    vector<pair<int,int>> frontier;
    for(int i = 0; i < game.size(); i++) {
        for(int j = 0; j < game[i].size(); j++) {
            if(game[i][j] >= 0) {
                for(pair<int,int> n : neighbors(i, j, game.size())) {
                    if(game[n.first][n.second] == -2) {
                        frontier.push_back(pair<int,int>(i,j));
                        break;
                    }
                }
            }
        }
    }
    return frontier;
}

vector<pair<int,int>> frontier_outer(vector<vector<int>>& game) {
    vector<pair<int,int>> frontier;
    for(int i = 0; i < game.size(); i++) {
        for(int j = 0; j < game[i].size(); j++) {
            if(game[i][j] == -2) {
                for(pair<int,int> n : neighbors(i, j, game.size())) {
                    if(game[n.first][n.second] != -2) {
                        frontier.push_back(pair<int,int>(i,j));
                        break;
                    }
                }
            }
        }
    }
    return frontier;
}

vector<vector<int>> eff_grid(vector<vector<int>>& grid) {
    vector<vector<int>> eg;
    eg = grid;
    for(int i = 0; i < eg.size(); i++) {
        for(int j = 0; j < eg[0].size(); j++) {
            if(grid[i][j] > 0) {
                for(pair<int,int> n : neighbors(i, j, eg.size())) {
                    if(grid[n.first][n.second] == -1)
                        eg[i][j]--;
                }
            }
        }
    }
    return eg;
}

vector<vector<int>> buildGrid(int size, double p, bool safe_origin) {
    vector<vector<int>> grid(size);
    for(int i = 0; i < size; i++)
        grid[i] = vector<int>(size, 0);
    
    int mines = (int)(p*size*size);

    vector<pair<int,int>> coordinates;
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            if(!safe_origin || (i > 1 && i < (size-1)) || (j > 1 && j < (size-1)))
                coordinates.push_back(pair<int,int>(i, j));
        }
    }

    vector<pair<int,int>> mine_coords;
    std::sample(coordinates.begin(), coordinates.end(), back_inserter(mine_coords), mines, mt19937{random_device{}()});

    for(pair<int,int> coord : mine_coords)
        grid[coord.first][coord.second] = -1;

    vector<vector<int>> new_grid(size);
    for(int i = 0; i < size; i++)
        new_grid[i] = vector<int>(size, 0);

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            if(grid[i][j] != -1) {
                for(pair<int,int> n : neighbors(i, j, size)) {
                    new_grid[i][j] -= grid[n.first][n.second];
                }
            }
            else
                new_grid[i][j] = -1;
        }
    }

    return new_grid;
}

void pretty_print(vector<vector<int>>& game) {
    map<int, char> ch;
    ch[-2] = '?'; ch[-1] = 'X'; ch[0] = ' ';
    for(int i = 1; i <= 8; i++)
        ch[i] = (char)(48+i);
    
    for(int i = 0; i < 6*game.size()+1; i++)
        printf("_");
    printf("\n");
    for(int i = 0; i < game.size(); i++) {
        printf("|");
        for(int j = 0; j < game[i].size(); j++)
            printf("     |");
        printf("\n|");
        for(int j = 0; j < game[i].size(); j++)
            printf("  %c  |", ch[game[i][j]]);
        printf("\n|");
        for(int j = 0; j < game[i].size(); j++)
            printf("_____|");
        printf("\n");
    }
}