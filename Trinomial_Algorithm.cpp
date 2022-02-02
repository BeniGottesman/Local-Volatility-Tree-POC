//std
#include <iostream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>

//QT-libraries
#include <QApplication>
#include <QVector>
#include <QFile>
#include <QDebug>

typedef QVector<float> v1fs;
typedef QVector<v1fs *> v2fs;

typedef QVector<float> v1f;
typedef QVector<v1f> v2f;

void plot_StoppingRegion ();
void simple_Pricing ();
void stoppingRegion(float n_abscisse = 100);
float LVDP (float y0);
float LVDP2 (float y0);
float LVDP_Binomial (float y0);
float LVDP_Trinomial (float y0, float T);
inline float option(float S, float K);
inline float put (float S, float K);
inline float call (float S, float K);
inline float p_BS (float X0, float S_0);
float p_AB (float X0, float S_0, int m = 100);
inline float p_sigmaJump (float X0, float S_0, float K);
//float exp_integral (float x0, float x, int n);
inline float A_BS (float x , float sqrt_h, float sigma_bar);
inline float A_AB (float x , float sqrt_h, float sigma_bar);
float Payoff (float t, float S, float K);
bool cmpf(float A_AB, float B, float epsilon = 0.005f);
inline float mu_CEV (float alpha_model, float beta_model, float z);
inline float sigma_CEV (float sigma_extern, float z);
inline float mu_AB (float z);
inline float sigma_AB (float z);
inline float Discount (float t);
float LVDP_Trinomial_Jump_Vol (float y0, float Maturity, float S1);//T = maturity
inline float A_sigma_jump (float x, float sqrt_h, float sigma_bar, float sigma_extern);
void table_simulation ();
inline float A_Dly (float x, float sqrt_h, float sigma_bar, float sigma);
inline float p_Dly (float X0, float S_0, float K);
float LVDP_Simple_Model (float S0, float K, float Penalty, float r, float sigma, float n, float Maturity);//T = maturity


float h;
/*const*/ int n = 1000;

//CEV Model Parameters
//float sigma_extern = 0.5;
float sigma_up = 1.0/30.0;
float sigma_down = 1.0/30.0;
float S_0 = 100.0;//==y
//float T = 0.5;
//float r = 0.1;
//float const B = 0.01;
//float const C = 200.0;
float const B = 0.10;//1.0=L;
float const C = 20.0;//20.0=U;
float beta_model = 2.0;
float alpha_model = 0.5;
float K = 100.0;
float Penalty = 12.0;//2.0
float sigma = 0.5;
//CEV Model Parameters

//AB Model Parameters
//float const y = 4.0;
float T = 2.0;
//float const T = 2.0;
float const r = 0.06;

//Other model
float const A1 = 10.0; float A2 = 10.0;
float const B1 = 2.0; float B2 = 2.0;
//Other model
//float const Strike = 4.0;
//AB Model Parameters

using namespace std;

int main()
{

//1 SIMPLE TEST//
//        float S0 = 90.0;
//        float vol = 0.4;
//        float Maturity = T;
//        float sqrt_h = sqrt(Maturity/n);
//        float sigma_bar = 0.5;//vol*C*sqrt(C)/30.0 + r*C*sqrt_h + 0.1f;
//        float price = LVDP_Simple_Model (S0, K, Penalty, r, sigma_bar, n, Maturity);
//        cout << "price = " << price << endl;
//1 SIMPLE TEST//

//STOPPING REGION//
    stoppingRegion(10);
//STOPPING REGION//

    return 0;
}


void stoppingRegion(float n_abscisse)
{
    float price_Trinomial, last_price_Trinomial;
    float tmp_payoff;
    float u, S0;
    float S0_max_ordinate = C-0.1;
    int m_S0 = 200;//nb of simulation on every ordinate axe
    float dS0 = abs(S0_max_ordinate-B)/m_S0;//ordinate division
    double elapsed_secs_Trinomial;
    float minValue, maxValue;
    float precision = 0.001f;
    float dt = T/n_abscisse;
    bool intoTheRegion;
    float S0_min, S0_max;
    float beginning_time = 0.0;
    string SR_down, SR_up;
    SR_down += "SR_down = fliplr([";
    SR_up += "SR_up = fliplr([";

    std::cout << "Stopping Region" <<endl;
    std::cout << "Parameters n = " << n << ", r = " << r << ", K = " << K << ", C = " << C
              << ", B = " << B << ", T = " << T << ", Penalty = " << Penalty << endl;
    std::cout << "Absciss = " << n_abscisse << ", Ordinate = " << m_S0 << endl;

    clock_t begin = clock();

    //Problem in final point ? at t=T value is buyer payoff (without penalty) check it ..
    SR_down += to_string(option(S0,K));
    SR_up += to_string(option(S0,K));
    SR_down += ",";
    SR_up += ",";

    for (int i = n_abscisse; i >= 0; i--)
    {
        //t = T - (float) (i*dt);//abscisse
        u = beginning_time + (float) (i*dt);//abscisse
        minValue = maxValue = 0.0;
        intoTheRegion = false;

        //reinitialization of the values
        price_Trinomial = last_price_Trinomial = 0.0;
        S0_min = S0_max = 0.0;
        //reinitialization of the values

        //cout << i << endl;

        for (int j = 0; j <= m_S0; j++)
        {
            S0 = 60.0f + j*dS0;

            price_Trinomial = LVDP_Simple_Model(S0, K, Penalty, r, sigma, n, T-u);//u=maturity
            tmp_payoff = option(S0,K)+Penalty;//max (K-S0, (float) 0.0) + Penalty;

            if (fabs(price_Trinomial-tmp_payoff) < precision)
            {
                if (!intoTheRegion)//We enter the Stopping Region
                {
                    S0_min = S0;
                    SR_down += to_string(S0);
                    SR_down += ",";
                    minValue = price_Trinomial;
                    intoTheRegion = true;
                    elapsed_secs_Trinomial = double(clock() - begin) / CLOCKS_PER_SEC;
                }
            }
            else if (intoTheRegion)//first upper point out of the region
            {
                maxValue = last_price_Trinomial;
                SR_up += to_string(S0_max);
                SR_up += ",";
                intoTheRegion = false;

                //!!! accelerate the simulation
                break;
                //!!!!
            }
            S0_max = S0;
            last_price_Trinomial = price_Trinomial;
        }
    }

    clock_t end = clock();
    elapsed_secs_Trinomial = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "n = " << n_abscisse << ", Time = " << elapsed_secs_Trinomial << " s, " <<
                 elapsed_secs_Trinomial/60 << " mn, " <<
                 elapsed_secs_Trinomial/(60*60) << " h" << endl;

    SR_down += "]);";
    SR_up += "]);";
    cout << "t = 0:"<<T<<"/"<<n_abscisse<<":"<< T << ";" << endl;
    cout << "hold on"<< endl;
    cout << SR_down << endl;
    cout << SR_up << endl;
    cout << "plot(t,SR_up,'DisplayName', 'Up frontier')" << endl;
    cout << "plot(t,SR_down,'DisplayName', 'Down frontier')" << endl;
    cout << "xlabel('Time');"<< endl;
    cout << "ylabel('Stopping Region');"<< endl;
    cout << "legend"<< endl;
}

void table_simulation ()
{
    double S1 = 100.0;
    std::cout << "Parameters $r = " << r << "$, $K = " << K << "$, $B = " << B
              << "$, $C = " << C << "$, $T = " << T << "$, $X0 = " << S_0 << "$, $\\delta = " << Penalty << endl;
    std::cout << ", $S_1 = " << S1 <<"$, $\\underline \\sigma = " << sigma_down << "$, $\\overline\\sigma = " << sigma_up
              << "$, $\\sigma= " << sigma << "$" << endl;

    double price;

    int tmp [4] = {2000, 4000, 8000, 10000};
    bool finish;
    double S0;
    for(int i = 0; i < 20; i++)
    {
        S0 = 80.0 + ((double) i)*5.0;
        cout << S0;
        for(int j = 0; j < 4; j++)
        {
            n  = tmp [j];
            price = LVDP_Trinomial_Jump_Vol(S0, T, S1);
            cout << " & " << price;
        }
        cout << "\\\\" << endl;
    }

}

void simple_Pricing ()
{
    clock_t begin = clock();
    float price_Trinomial = LVDP_Trinomial_Jump_Vol(S_0, T, 100.0);//LVDP_Trinomial (y, T);
    clock_t end = clock();

    double elapsed_secs_Trinomial = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Parameters r = " << r << ", K = " << K << ", C = " << C
              << ", B = " << B << ", T = " << T << ", X0 = " << S_0 << ", Penalty = " << Penalty << endl;
    std::cout << "Trinomial Price = " << price_Trinomial << ", n = " << n << ", Time = " << elapsed_secs_Trinomial << " s" << endl;
    std::cout << "sigma down = " << sigma_down << ", sigma up= " << sigma_up << ", sigma= " << sigma << endl;
}

//1st main algorithm
float LVDP_Simple_Model (float S0, float K, float Penalty, float r, float sigma, float n, float Maturity)//T = maturity
{
    h = Maturity/n;//time unity
    float tmp_f, tmp_g, tmp_dp;
    float sqrt_h = sqrt(h);
    float sum_xi_kj;
    float q_up, q, q_down;//Probabilities

    float tmp_t = (float (n)) * h;
    float Discount_t = exp (-r*tmp_t);
    float sigma2 = sigma*sigma;
    float sig2, tmp1, tmp2, tmp;

    int nb_of_leafs = 2*n + 1;
    int nb_of_nodes;// = 0.0;
    v1f * J_t_plus_1 = new v1f (nb_of_leafs);
    v1f * J_t;

    float S, sig;
    //INITIALIZATION
    //#pragma omp parallel for
    for (int j = 0; j < nb_of_leafs; j++)
    {
        sum_xi_kj = float (n-j);
        S = S0*exp(sigma*sqrt_h*(sum_xi_kj));
        (*J_t_plus_1) [j] = max(0.0f, S - Discount_t*K);//Payoff (tmp_t, z, K);
    }

    for (int time = n-1; time >= 0; time--)
    {
        tmp_t = (float (time)) * h;
        Discount_t = exp (-r*tmp_t);//because exponential is very slow
        nb_of_nodes = 2*time + 1;
        J_t = new v1f (nb_of_nodes);

        for (int j = 0; j < nb_of_nodes; j++)
        {
            S = S0*exp(sigma*sqrt_h*(j-nb_of_nodes));
            sig = max(0.05, min(0.5, sqrt(S)/30.0));
            sig2 = sig*sig;

            q_up = 0.5*( sig2/(sigma2) - sig2*sqrt_h/(2.0*sigma) );
            q_down = sig2/(sigma2) - q_up;
            q = 1.0-q_down-q_up;

            tmp = q_up*(*J_t_plus_1)[j+2]
                    + q*(*J_t_plus_1)[j+1]
                    + q_down*(*J_t_plus_1)[j];
            float tmp_option = max(0.0f, S - Discount_t*K);
            tmp1 = max (tmp_option, tmp);
            tmp2 = max (0.0f, tmp_option + Discount_t*Penalty);
            (*J_t) [j]= min (tmp1, tmp2);
        }
        delete J_t_plus_1;
        J_t_plus_1 = J_t;
    }

    float price = (*J_t)[0];
    delete J_t, J_t_plus_1;

    return price;
}

//2nd Algorithm with Jumps
//S1 is for the volatility
float LVDP_Trinomial_Jump_Vol (float y0, float Maturity, float S1)//T = maturity
{
    h = Maturity/n;//time unity
    float tmp_f, tmp_g, tmp_dp;
    float sqrt_h = sqrt(h);
    float sum_xi_kj;
    float q_m1, q_0, q_p1;//Probabilities

    //WORK BETTER !!!!
    float BC = max (abs(B), abs(C));
    float sigma = 1.0/30.0*sqrt(BC);
    //jump volatility
    //float sigma = max (sigma_up, sigma_down);
    float sigma_bar= sigma*BC + r*BC*sqrt_h + .1f;
    //WORK BETTER !!!!
    
    float tmp_sigma, tmp_mu;
    float X0 = y0;

    int nb_of_leafs = 2*n + 1;
    int nb_of_nodes;// = 0.0;
    v1f * J_t_plus_1 = new v1f (nb_of_leafs);
    v1f * J_t;


    float tmp_t = (float (n)) * h;
    float exponent;
    float Discount_t = exp (-r*tmp_t);//We calculate once because exponential is very slow
    
    float z;
    float bn = min ( (float) n, (X0-B)/(sigma_bar*sqrt_h) );
    float cn = min ( (float) n, (C-X0)/(sigma_bar*sqrt_h) );

    //INITIALIZATION
    //#pragma omp parallel for
    for (int j = 0; j < nb_of_leafs; j++)
    {
        sum_xi_kj = float (n-j);
        z = X0 + sqrt_h*sigma_bar*(max(-bn, min(cn, sum_xi_kj)));
        (*J_t_plus_1) [j] = Discount_t*option(z,K);//Payoff (tmp_t, z, K);
    }

    float p_msh, p_psh, tmp_A, p_mA, p_pA, p_X0;
    for (int time = n-1; time >= 0; time--)
    {
        tmp_t = (float (time)) * h;
        Discount_t = exp (-r*tmp_t);//because exponential is very slow
        nb_of_nodes = 2*time + 1;
        J_t = new v1f (nb_of_nodes);

        for (int j = 0; j < nb_of_nodes; j++)
        {
            sum_xi_kj = float (time-j);
            z = X0 + sqrt_h*sigma_bar*(max(-bn, min(cn, sum_xi_kj)));
            tmp_f = Discount_t*option(z,K);//max ( (K - z), (float) 0.0 );// * Payoff (t, z, K);
            tmp_g = tmp_f + Discount_t*Penalty;
            {

                //IF JUMP VOL : sigma//
                //                if(z > S1)
                //                    tmp_sigma = sigma_down;
                //                else //if (z < S1)
                //                    tmp_sigma = sigma_up;
                //IF JUMP VOL : sigma//

                //Our model
                tmp_sigma = z*sqrt(z)/30.0;
                tmp_mu = z*r;
                //Our model

                exponent = -2.0f*tmp_mu/tmp_sigma;

                //Continuous Model uncomment//
                //                tmp_A = A_Dly(z, sqrt_h, sigma_bar, tmp_sigma);//z*z/sigma_bar*sqrt_h;
                //                p_msh = p_Dly(X0, z - sigma_bar*sqrt_h, exponent);
                //                p_psh = p_Dly(X0, z + sigma_bar*sqrt_h, exponent);
                //                p_mA = p_Dly(X0, z-tmp_A, exponent);
                //                p_pA = p_Dly(X0, z+tmp_A, exponent);
                //                p_X0 = p_Dly(X0, z, exponent);
                //Continuous Model uncomment//
                
                //Jump Model uncomment//
                //                tmp_A = A_sigma_jump(z, sqrt_h, sigma_bar, tmp_sigma);//z*z/sigma_bar*sqrt_h;
                //                p_msh = p_sigmaJump(X0, z - sigma_bar*sqrt_h, exponent);
                //                p_psh = p_sigmaJump(X0, z + sigma_bar*sqrt_h, exponent);
                //                p_mA = p_sigmaJump(X0, z-tmp_A, exponent);
                //                p_pA = p_sigmaJump(X0, z+tmp_A, exponent);
                //                p_X0 = p_sigmaJump(X0, z, exponent);
                //Jump Model uncomment//
                
                //Probabilities
                q_p1 = (p_X0-p_mA)*(p_pA-p_X0) / ( (p_pA-p_mA)*(p_psh-p_X0) );
                q_m1 = (p_pA-p_X0)*(p_X0-p_mA) / ( (p_pA-p_mA)*(p_X0-p_msh) );
                q_0 = 1.0f - q_m1 - q_p1;
                //Probabilities

                tmp_dp = q_p1*(*J_t_plus_1)[j] +
                        q_0*(*J_t_plus_1)[j+1] +
                        q_m1*(*J_t_plus_1)[j+2];
                (*J_t) [j] = min (tmp_g, max (tmp_f, tmp_dp));//Game Option
                //(*J_t) [j] = max (tmp_f, min (tmp_g, tmp_dp));//Game Option
                //(*J_t) [j] = max (tmp_f, tmp_dp);//American Option
            }
        }
        delete J_t_plus_1;
        J_t_plus_1 = J_t;
    }

    float price = (*J_t)[0];
    delete J_t, J_t_plus_1;

    return price;
}

inline float A_Dly (float x, float sqrt_h, float sigma_bar, float sigma)
{
    return x*x*sigma*sigma/sigma_bar * sqrt_h;
}

inline float p_Dly (float X0, float y, float K)
{
    return 0.0;//x*x*sigma*sigma/sigma_bar * sqrt_h;
}

inline float A_sigma_jump (float x, float sqrt_h, float sigma_bar, float sigma)
{
    return x*x*sigma*sigma/sigma_bar * sqrt_h;
}

inline float p_sigmaJump (float X0, float y, float K)
{
    return pow(y,K+1)/pow(X0,K)-X0;
}

inline float p_BS (float X0, float y)
{
    return X0-X0*X0/y;
}

inline float A_BS (float x, float sqrt_h, float sigma_bar)
{
    return x*x/sigma_bar * sqrt_h;
}

float p_AB (float X0, float y, int m)
{
    float I_value = 0.0;
    float exp_I = 0.0;
    float sigma_square;
    float tmp_exp_Ivalue = 0.0;
    float x;
    float dz = (y-X0)/((float) m);

    for (int i = 0; i < m; i++)
    {
        x = X0 + i*dz;
        sigma_square = sigma_AB(x)*sigma_AB(x);
        tmp_exp_Ivalue += (mu_AB(x)/sigma_square)*dz;
        exp_I += exp (-2.0 * tmp_exp_Ivalue);//exp is very slow
        I_value += exp_I*dz;
    }

    return I_value;
}


inline float A_AB (float x , float sqrt_h, float sigma_bar)
{
    return sigma_AB(x)*sigma_AB(x)/(sigma_bar)*sqrt_h;
}

//CEV model
inline float mu_CEV (float alpha_model, float beta_model, float z)
{
    return beta_model - alpha_model*z;
}

inline float sigma_CEV (float sigma, float z)
{
    return sigma*sqrt (z);
}
//CEV model

//AB model
inline float mu_AB (float z)
{
    return max(B1, min (A1, z));
}

inline float sigma_AB (float z)
{
    return max(B2, min (A2, z));
}
//AB model

inline float Payoff (float t, float z, float K)
{
    if (K - z > 0.0)
        return /*exp(-r*(t))**/(K - z);// but exponential is very slow to compute each time
    return 0.0;
}

inline float Discount (float t)
{
    return exp(-r*(t));
}

bool cmpf(float A, float B, float epsilon)
{
    return (abs(A - B) < epsilon);
}

inline float option (float S, float K)
{
    return call (S,K);
}

inline float put (float S, float K)
{
    return max (0.0f, K-S);
}

inline float call (float S, float K)
{
    return max (0.0f, S-K);
}
