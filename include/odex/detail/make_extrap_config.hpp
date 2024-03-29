
#ifndef ODEX_DETAIL_MAKE_EXTRAP_CONFIG_HPP
#define ODEX_DETAIL_MAKE_EXTRAP_CONFIG_HPP

#include <algorithm>
#include <cassert>
#include <vector>
#include <tuple>

namespace odex {
namespace detail {

template <class T>
inline auto make_extrap_config(std::size_t order, std::size_t num_cores)
{
    float isbn = 0.0f;
    std::vector<std::size_t> step_counts;
    std::vector<long double> weights;

    if (order == 8)
    {
        if (num_cores == 3)
        {
            isbn = 0.5799f;
            step_counts = {
                2, 16, 18, 20
              };
            weights = {
               -1.0L/498960.0L,
                65536.0L/9639.0L,
               -531441.0L/25840.0L,
                250000.0L/16929.0L
              };
        }
        else if (num_cores == 6)
        {
            isbn = 0.7675f;
            step_counts = {
                2, 4, 6, 10, 8, 12, 14, 16, 18, 20, 22
              };
            weights = {
               -32952289146985386285870523118228405533963.0L/8936455970950449255004500793755553651752960000.0L,
                577598451788090848795408620332945866052063.0L/7941577083559481271537202853825736155366400000.0L,
                85250432905463981456535914913119571901637.0L/122585129917015764814876554098155742822400000.0L,
                1677712357266484804784340039643670407130779.0L/200176613749290063312100817780124401799266304.0L,
                2165.0L/767488.0L,
                13805.0L/611712.0L,
                4553.0L/72080.0L,
                14503.0L/66520.0L,
                27058.0L/7627.0L,
               -86504.0L/5761.0L,
                40916.0L/3367.0L
              };
        }
        else if (num_cores == 8)
        {            
            isbn = 0.8176f;
            step_counts = {
                2, 26, 28, 30, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24
              };
            weights = {
               -298857882660976887631476729981565763568862608650111.0L/418309165211319520505929581345807932301941968522444800000.0L,
                54841752514603990885070634946141665271319680054382001869.0L/7796886807193233666234782137510621223379391720980480000.0L,
               -6653387365673258947809103108875129803987861502988566763111.0L/258933840128714385112278118363748576398397951218483200000.0L,
                54824130826438857272172198804804549641875992497090297.0L/2867295913488162504174863944731843756709312047611904.0L,
                6833.0L/476577792.0L,
                10847.0L/91078656.0L,
                15235.0L/34643968.0L,
                383.0L/321152.0L,
                543.0L/198784.0L,
                9947.0L/1741056.0L,
                6243.0L/543104.0L,
                6875.0L/296192.0L,
                1401.0L/28496.0L,
                17713.0L/152688.0L,
                6375.0L/19264.0L
              };
        }
    }
    else if (order == 12)
    {
        if (num_cores == 4)
        {
            isbn = 0.4515f;
            step_counts = {
                2, 8, 12, 14, 16, 20
              };
            weights = {
               -1.0L/157172400.0L,
                4096.0L/155925.0L,
               -59049.0L/15925.0L,
                282475249.0L/15752880.0L,
               -4194304.0L/178605.0L,
                9765625.0L/954261.0L
              };
 
        }
        else if (num_cores == 8)
        {
            isbn = 0.7116f;
            step_counts = {
                2, 8, 10, 16, 24, 26, 4, 6, 12, 14, 18, 20, 22, 28, 30
              };
            weights = {
               -1703338201142081344537976944145527211643949659234240721389419.0L/23648864513368626787371236562816879339803777703368508907192320000000000.0L,
                28566269141029842679611128435317644430416456404930682840133.0L/1235974431889711160110009091223554591673172898357667840000000000.0L,
                1661823701099033749417849761031734684833334503871915993221173.0L/16039458446054067082385395561773359826518343984165155963236515840.0L,
                297002124618857676974925717053765105019453996390390125558609.0L/160179791893258872271743935365682835875612617239142743750000000.0L,
               -5460019744535790351900106662607930219497507008045052153266932061.0L/109934733569605065449891190520737372992080772483328000000000000.0L,
                4518788471550054059819510090434891452487764271627191207619322033987247547.0L/24806501237799258867871926464493230076717249339197736615936000000000000.0L,
                235.0L/21030240256.0L,
                4147.0L/1612709888.0L,
                11521.0L/39731200.0L,
                2375.0L/3528704.0L,
                6435.0L/708736.0L,
                1291.0L/15780.0L,
                11311.0L/4672.0L,
               -180864.0L/751.0L,
                222080.0L/2079.0L
              };
        }            
    }
    else if (order == 16)
    {
        if (num_cores == 5)
        {
            isbn = 0.4162f;
            step_counts = {
                2, 8, 10, 12, 14, 16, 18, 22
              };
            weights = {
               -1.0L/365783040000.0L,
                4194304.0L/456080625.0L,
               -6103515625.0L/11955879936.0L,
                544195584.0L/74449375.0L,
               -678223072849.0L/17079828480.0L,
                68719476736.0L/749962395.0L,
               -2541865828329.0L/31682560000.0L,
                379749833583241.0L/16878274560000.0L
              };
        }
    }
    assert(step_counts.size() != 0 && "No matching extrapolation scheme for order and number of cores!");

    std::vector<T> retweights(weights.size());
    std::transform(weights.begin(), weights.end(), retweights.begin(),
        [](long double d){ return static_cast<T>(d); });
    return std::make_tuple(isbn, step_counts, retweights);
}

} // namespace detail
} // namespace odex

#endif // ODEX_DETAIL_MAKE_EXTRAP_CONFIG_HPP

