/**
*   @file lozi.cpp
*   @name Mapas de Lozi
*   @author Julio Lugo
*   @mail jmanuellugo96@gmail.com
*   @date 16/05/2017
*   @brief Globaly coupled Lozi map 
*/

# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>
# include <ctime>
# include <cstdlib>
# include <cstdio>
# include <random>
# include <vector>
# include <thread>

using namespace std;

# define alpha 1.4 ///Alpha parameter (fixed)
# define b 0.5 ///B fixed parameter


///How to use the program

void usage()
{
  cout << "Execute as follows: \n\n";
  cout << "\t[./your_chosen_name.exe],\n\t[epsilon value]\n\t[beta value],\n\t[mu value],\n\t[population i size],\n\t[population j size],\n\t[Total iterations],\n\t[Threshold]" << endl;
}




///f(x,y) function definition
inline double function_xy(const double &x, const double &y)
{
  return 1 - alpha * fabs(x) + y;
}

inline double function_x(const double &x)
{
  return  (1 - pow(b, (1 - x)*x))/(1 - sqrt(sqrt(b)));
}

///x value function
inline double x_function(const double &epsilon, const double &mu,
                        const double &f, const double &h_i, const double &h_j)
{ 
  return (1 - epsilon)*f + epsilon * h_j + mu * h_i;
}


///y value function
inline double y_function(const double &x, double &beta)
{
  return beta * x;
}






///Routine for initial conditions setting
void set_initial_conditions(double **x_axis, double **y_axis,
                            double **f_var, const int * num_nodes, char &s)
{
    random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

    if(s != '4')
    {
      uniform_real_distribution<> dis(-1, 1);

      for(int i = 0; i < 2; i++)
      for(int j = 0; j < num_nodes[i]; j++)
        {
          x_axis[i][j] = dis(gen);
          y_axis[i][j] = dis(gen);
        }
    }
    else
    {
      uniform_real_distribution<> dis(0, 1);

      for(int i = 0; i < 2; i++)
      for(int j = 0; j < num_nodes[i]; j++)
        {
          x_axis[i][j] = dis(gen);
          y_axis[i][j] = dis(gen);
        }
    }



    

}






/// Copy values from x_axis to m_axis, it is used in threshold function

inline void copy_values(double **x_axis, double **m_axis, const int *num_nodes)
{
  for (int i = 0; i < 2; i++)
    for(int j = 0; j < num_nodes[i]; j++)
      m_axis[i][j] = x_axis[i][j];
}



///Set min value allows to sort all pops

void set_min_value(double **x_axis, double **m_axis,
                    const int *num_nodes)
{
  int row = 0, col = 0;
    

  for(int i = 0; i < 2; i++)
    for (int j = 0; j < num_nodes[i]; j++)
      {
        double min = 2.0;

        for (int k = 0; k < 2; ++k)
          for(int l = 0; l < num_nodes[k]; l++)
            if(m_axis[k][l] < min)
              {
                min = m_axis[k][l];
                row = k;
                col = l;
              }

        m_axis[row][col] = 2.0;
        x_axis[i][j] = min;

      }
}


void check_clusters(double **x_axis, const int *num_nodes)
{
  vector < tuple<double, int, int> > values[2];
  tuple<double, int, int> aux;
  pair <double, int> special_a(0,1), special_b(0,0);
  bool found = false;

  for(int i = 0; i < 2; i++)
    for (int j = 0; j < num_nodes[i]; j++)
    {
      found = false;

      for(auto &p : values[i])
        if(x_axis[i][j] == get<0>(p))
        {
          get<1>(p)++;
          found = true;
          break;
        }

      if(not found)
        values[i].push_back(make_tuple(x_axis[i][j], 1, j));
    }
/*
   for(int i = 0; i < 2; i++)
   {
      

      cout << "Pop " << i << endl;

      for(auto &p : values[i])
        {
          cout << get<0>(p) << " cant: " << get<1>(p) << " pos: " << get<2>(p) << endl;
        }

     cout << "\n";
   }*/

    if(values[0].size() == 1 and values[1].size() == 1)
      return;

    found = false;

  if(fabs(get<0>(values[0].back()) - get<0>(values[1].front()) < 0.05))
    {
     //   cout << "Condición especial!" << endl;

        for(int j = num_nodes[0] - 1; j > 0; j--)
        {
          double result = fabs(x_axis[0][j] - x_axis[0][j - 1]);
          double min = 0.05;



          if(result < min)
          {

    //      cout << x_axis[0][num_nodes[0] - 1] << " " << x_axis[0][j-1] << endl;
            special_a.second++;
          }
          else
          {
            
            special_a.first = x_axis[0][j - 1];
            break;
          }
        }


     //  cout << "\n\n";
        for(int j = 0; j < num_nodes[1] - 1; j++)
        {
          double result = fabs(x_axis[1][j] - x_axis[1][j + 1]);
          double min = 0.05;

          if(result < min)
          {

     //     cout << "diff: " << result << endl;
     //     cout << "" << x_axis[0][j] << " " << x_axis[1][j + 1] << endl;
            special_b.second++;
          }
          else
          {
            special_b.first = x_axis[1][j + 1];
            break;
          }
        }

    //    cout << "pop a: " << special_a.second << " pop b: " << special_b.second << endl;

        if(special_a.second > special_b.second)
        {
          if(get<1>(values[1][1]) > 1)
            for(int j = 0; j <= special_b.second; j++)   
              x_axis[1][j] = get<0>(values[1][1]);
          else
            for(int j = 0; j <= special_b.second; j++)
            {
              random_device rd;  //Will be used to obtain a seed for the random number engine
              mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
              uniform_real_distribution<> dis(special_b.first - 0.0001, special_b.first + 0.0001);

              x_axis[1][j] = dis(gen);

            }
        }
        else
        {
          if(get<1>(values[0][values[0].size()-1]) > 1)
            for(int j = num_nodes[0] - 1; j > (num_nodes[0] - 1 - special_a.second); j--)  
              x_axis[0][j] = get<0>(values[0][values[0].size() - 2]);
          else
            for(int j = num_nodes[0] - 1; j > (num_nodes[0] - 1 - special_a.second); j--)
            {
              random_device rd;  //Will be used to obtain a seed for the random number engine
              mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
              uniform_real_distribution<> dis(special_a.first - 0.0001, special_a.first + 0.0001);

              x_axis[0][j] = dis(gen);
            }
        }
    }

      values[0].clear();

      values[1].clear();


  for(int i = 0; i < 2; i++)
    for (int j = 0; j < num_nodes[i]; j++)
    {
      found = false;

      for(auto &p : values[i])
        if(x_axis[i][j] == get<0>(p))
        {
          get<1>(p)++;
          found = true;
          break;
        }

      if(not found)
        values[i].push_back(make_tuple(x_axis[i][j], 1, j));
    }
/*
   for(int i = 0; i < 2; i++)
   {
      

      cout << "AFTER OPERATION " << i << endl;

      for(auto &p : values[i])
        {
          cout << get<0>(p) << " cant: " << get<1>(p) << " pos: " << get<2>(p) << endl;
       }

      cout << "\n";
   }

*/


}

///Threshold function: it allows to sort the population


void threshold_function(double **x_axis, double **m_axis, const int *num_nodes)
{

  copy_values(x_axis, m_axis, num_nodes);

  set_min_value(x_axis, m_axis, num_nodes);

  check_clusters(x_axis, num_nodes);
}





///Print all both populations in one file

inline void plot_print(int t, ofstream &out, double **x_axis, 
                        const int *num_nodes)
{
  for(int i = 0; i <num_nodes[0]; i++)
    out <<  fixed << setprecision(5) << x_axis[0][i] << "\t";

  for(int i = 0; i < num_nodes[1]; i++)
    out <<  fixed << setprecision(5) << x_axis[1][i] << "\t";

  out << "\n";

}

inline int theta(double result)
{
  if(result < 0)
      return 0;
  else
      return 1;
}

pair < double, double > p_function(double **x_axis, const int * num_nodes)
{
  double min_delta = 0.000005, distance = 0, productory = 1, summatory = 0;
  pair <double, double> ret_val;


  for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < num_nodes[i]; ++j)
        {
          productory = 1;

          for(int k = 0; k < num_nodes[i]; k++)
          {
            if(j == k)
              continue;

            distance = fabs(x_axis[i][j] - x_axis[i][k]);

            productory *= theta(distance - min_delta);

          }

          summatory += productory;
        }

      if(i == 0)
        ret_val.first = 1 - summatory/num_nodes[i];
      else
        ret_val.second = 1 - summatory/num_nodes[i];

      summatory = 0;
    }


  return ret_val;
}



///Simulation main module. It returns a tuple which contains delta value and sigma from pop i and pop j. 

tuple <double, double *> simulate(double &epsilon, 
                                  double &beta, 
                                  double &mu, 
                                  int &pop_i, 
                                  int &pop_j, 
                                  const unsigned int &total_it, 
                                  unsigned int &threshold, 
                                  char &s)
{


  const int num_nodes[2] = {pop_i, pop_j}; ///Vector used to define population sizes
  const unsigned int tau = total_it - 1000;

  if(epsilon == 0.0)
    epsilon = 0.05;
  
  ofstream out1("pop.dat"); ///Archivo final
  
  double N = num_nodes[0] + num_nodes[1]; /// N = Ni + Nij

  double  ret_val_delta = 0, ret_val_sigma[2] = {0}, **y_axis = new double*[2], **x_axis = new double*[2],** f_var = new double*[2], **m_axis = new double*[2];

  for(int i = 0; i < 2; i++)
  {
    x_axis[i] = new double[num_nodes[i]];
    y_axis[i] = new double[num_nodes[i]];
    m_axis[i] = new double[num_nodes[i]];    
    f_var[i] = new double[num_nodes[i]];
  }

  set_initial_conditions(x_axis, y_axis, f_var, num_nodes, s);



  for(unsigned int t = 0; t < total_it; t++) ///Transient
  {  


    double h[2] = {0}; ///For each t, medium field is initialized


    for(int j = 0; j < 2; j++)
    {
      for(int k = 0; k < num_nodes[j]; k++)
      {      
        

       if(s == '4')
         f_var[j][k] = function_x(x_axis[j][k]);
      else
          f_var[j][k] = function_xy(x_axis[j][k], y_axis[j][k]);

        h[j] += f_var[j][k];
        y_axis[j][k] = y_function(x_axis[j][k], beta);

      }

      h[j] /= (N + 0.0);    ///Get medium field of time t

    }

    ///Set X values
    
    for(int j = 0, l = 1; j < 2; j++, l--)
      for(int k = 0; k < num_nodes[j]; k++)
        x_axis[j][k] = x_function(epsilon, mu, f_var[j][k], h[j], h[l]);
        

    if(t > threshold and t < threshold + 2)
      threshold_function(x_axis, m_axis, num_nodes);


    if(t > tau) ///If tau iteration is reached, then start printing
    {
       if(s != '2') ///Condition if you're diagramming or plotting
        plot_print(t + 1, out1, x_axis, num_nodes);

      ///Sigma and delta calculations
        double x_mean[2] = {0}, delta = 0, sigma[2] = {0};

      //  cout << "Tiempo: " << t << endl;

        for(int j = 0; j < 2; j++)
        {

            for(int k = 0; k < num_nodes[j]; k++)
            {
     //         cout << setprecision(5) << x_axis[j][k] << " ";
              x_mean[j] += x_axis[j][k];
            }
     //       cout << "\n\n";
            x_mean[j] /= num_nodes[j];
                  
        }



     //   cout << "Valor medio de Xi: " << setprecision(5) << x_mean[0] << endl;
       // cout << "Valor medio de Xj: " << setprecision(5) << x_mean[1] << endl;

        for(int j = 0 ; j < 2; j++)
        {
          for(int k = 0; k < num_nodes[j]; k++)
            sigma[j] += (x_mean[j] - x_axis[j][k])*(x_mean[j] - x_axis[j][k]);

            sigma[j] /= num_nodes[j];

            sigma[j] = sqrt(sigma[j]);

            ret_val_sigma[j] += sigma[j];
        }
        
   //     cout << "Sigma: " << setprecision(5) << sigma[0] << endl;
     //   cout << "Sigma: " << setprecision(5) << sigma[1] << endl;

        delta = fabs(x_mean[0] - x_mean[1]);

        ret_val_delta += delta;

    }


  }

    ret_val_delta /= (total_it - tau);

    if(ret_val_delta < 0.0001)
      ret_val_delta = 0;

    for(int j = 0; j < 2; j++)
    {
     ret_val_sigma[j] /= (total_it - tau);

     if(ret_val_sigma[j] < 0.025)
      ret_val_sigma[j] = 0;
    }

    if(s == '3')
    {
      auto ret_val_p = p_function(x_axis, num_nodes);
      ret_val_sigma[0] = ret_val_p.first;
      ret_val_sigma[1] = ret_val_p.second;
    }


  out1.close();


    for(int i = 0; i < 2; i++)
    {
      delete[] x_axis[i];
      delete[] y_axis[i];
      delete[] m_axis[i];
      delete[] f_var[i];
    }

  return make_tuple(ret_val_delta, ret_val_sigma);
}


///Initial menu displayed 

inline char menu()
{
  char ret_val;
  system("clear");
  cout << "\t\tGLOBALLY COUPLED LOZI MAP" << endl;

  cout << "\t\t-------------------------\n\n" << endl;
  cout << "\tPlease choose an option.\n\n" << endl;
  cout << "\t1 .- Plot with delta and sigma (i, j)" << endl;
  cout << "\t2 .- Get phase diagram at range [mu: 0] [beta: 0.4]" << endl;
  cout << "\t3 .- Plot with delta and p (i, j)" << endl;
  cout << "\t4 .- Plot with new equation (using sigma (i, j) )\n\n" << endl;
  cout << "PRESS CTRL + C to finish the program" << endl;
  cout << "\n\nYour selection: ";
  cin >> ret_val;

  return ret_val;
}


void plot_phase_diagram(double &epsilon, 
                        double &beta, 
                        double &mu, 
                        int &pop_i, 
                        int &pop_j, 
                        const unsigned int &total_it, 
                        unsigned int &threshold, 
                        int &rows,
                        int &cols)
{
  ofstream m("plot.dat"), d("coord.dat");

  vector <tuple <double, double*> > data; ///Here will be saved current delta and sigma values

  char s = 'n';
  double ep_prop = (0.4)/(rows - 1 + 0.0);
  double beta_prop = (0.8)/(cols - 1 + 0.0);


  double new_mu = 0, new_beta = 0, act_ep = 0;

  cout << "\n\n\tPaso de mu: " << ep_prop << "\tPaso de beta: " << beta_prop << endl;

  for(int i = 0; new_mu >= 0.0 and i < rows; i++)
    {

      new_mu = mu - i*ep_prop;
      act_ep = mu - i*ep_prop;

      if(new_mu < 0.05)
      {      
        new_mu = 0.05;
      }

      if(new_mu < 0.0)
        break;

      for(int j = 0; new_beta <= 0.40 and j < cols; j++)
        {
          new_beta = beta + j * beta_prop;
          double delta = 0, sigma_i = 0, sigma_j = 0;
          vector < tuple<double, double*> > aux;


          if(fabs(new_beta) < 0.000001)
            new_beta = 0.0;

          if(new_beta > 0.40)
            break;

          cout << "\t\t\tMu: " << act_ep << "\tBeta: " << new_beta << endl;
          cout << "\t\t--------------------------------------------------" << endl;



          for(int i = 0; i < 50; i++)
          {
            


            aux.push_back(simulate(epsilon, new_beta, new_mu, pop_i, pop_j, total_it, threshold, s));

            delta += get<0>(aux[i]);
            sigma_i += get<1>(aux[i])[0];
            sigma_j += get<1>(aux[i])[1];

          /*  cout << "Simu Nº: " << setw(4) << left << i + 1;
            cout <<  "Delta: "<< setw(12) << left << fixed << setprecision(5) << get<0>(aux[i]);
            cout << "Sigma i: " << setw(12) << left << fixed << setprecision(5) << get<1>(aux[i])[0];
            cout << "Sigma j: "<< setw(12) << left << fixed << setprecision(5) << get<1>(aux[i])[1] << endl;*/
              
          }

          delta /= 50.0;
          sigma_i /= 50.0;
          sigma_j /= 50.0;

          cout << "--------------------------------------------------" << endl;
          cout << "\n\nCoordinate: (" << new_beta << ", " << act_ep << "):\n"<< "\t\tDelta: "<< delta << "\tMean Sigma i: " << sigma_i << "\tMean Sigma j: " << sigma_j << "\n\nFINAL MATRIX VALUE: ";

            d << new_beta << "\t" << act_ep << "\t";

            if(delta == 0 and sigma_i == 0 and sigma_j == 0)
            {
              m << "0\t";
              d << "0" << endl;
              cout << "0\n\n" << endl;
            }
            else if(delta != 0 and sigma_i != 0 and sigma_j != 0)
            {
              m << "1\t";
              d << "1" << endl;
              cout << "1\n\n" << endl;
            }
            else if(delta != 0 and (sigma_i == 0 and sigma_j == 0))         
            {
              m << "2\t";
              d << "2" << endl;
              cout << "2\n\n" << endl;
            }
            else
            {
              m << "3\t";
              d << "3" << endl;
              cout << "3\n\n" << endl;
            }

            
        }

        m << "\n";

        new_beta = beta;
    }


}





int main(int argc, char *argv[])
{
  if(argc != 8)
  {
    usage();
    return 0;
  }

  double epsilon = atof(argv[1]), beta = atof(argv[2]), mu = atof(argv[3]); ///Lectura del parámetro de acoplamiento

  int pop_i = atoi(argv[4]), pop_j = atoi(argv[5]);

  const unsigned int total_it = atoi(argv[6]);

  unsigned int threshold = atoi(argv[7]); 

  if(epsilon < 0 or epsilon > 1)
  {
    cerr << "Invalid value for epsilon" << endl;
    return 0;
  }

  if(fabs(beta) > 0.40)
  {
    cerr << "Invalid value for beta" << endl;
    return 0;
  }

  if(fabs(mu) > 0.40)
  {
    cerr << "Invalid value for mu" << endl;
    return 0;
  }

  if(pop_j <= 0 or pop_i <= 0)
  {
    cerr << "One population has an invalid size" << endl;
    return 0;
  }

  if(threshold > total_it)
  {
    cerr << "Threshold must be lower than total iterations" << endl;
    return 0;
  }

  char s = menu();
  tuple <double, double*> p;
  int rows = 0, cols = 0;

  switch(s)
  {
    case '1':
            p = simulate(epsilon, beta, mu, pop_i, pop_j, total_it, threshold, s);
            cout << "Delta: " << get<0>(p) << "\t\tSigma i: " << get<1>(p)[0] << "\tSigma j: " << get<1>(p)[1] << endl;
            system("./script.sh");
            break;
    case '2':
            system("clear");
            cout << "\n\n\n\t\tIntroduzca cantidad de filas de su matriz: ";
            cin >> rows;
            cout << "\t\tIntroduzca cantidad de columnas de su matriz: ";
            cin >> cols;
            cin.ignore();
            plot_phase_diagram(epsilon, beta, mu, pop_i, pop_j, total_it, threshold, rows, cols);
            break;
    case '3':
            p = simulate(epsilon, beta, mu, pop_i, pop_j, total_it, threshold, s);
            cout << "Delta: " << get<0>(p) << "\t\tP i: " << get<1>(p)[0] << "\tP j: " << get<1>(p)[1] << endl;
            system("./script.sh");
            break;
    case '4':
            p = simulate(epsilon, beta, mu, pop_i, pop_j, total_it, threshold, s);
            cout << "Delta: " << get<0>(p) << "\t\tSigma i: " << get<1>(p)[0] << "\tSigma j: " << get<1>(p)[1] << endl;
            system("./script.sh");
            break;
    default: break;

  }

  return 0;
}