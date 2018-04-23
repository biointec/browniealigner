#include <gtest/gtest.h>
#include <cstdio>
#include "util.h"

using namespace std;

TEST(Util, MixtureModelTest)
{
      /*  map<unsigned int, double> data;

        // read the data
        ifstream ifs("spectrum.txt");
        while (true) {
                double x, y;
                ifs >> x >> y;

                if (!ifs.good())
                        break;

                data[(unsigned int)x] = y;
        }

        ifs.close();

        // initialize the variables
        vector<double> mu(3), var(3), MC(3);

        mu[0] = 2;
        var[0] = 0;
        MC[0] = 1;

        mu[1] = 145;
        var[1] = 450;
        MC[1] = 1;

        mu[2] = 681;
        var[2] = 900;
        MC[2] = 1;

        Util::binomialMixtureEM(data, mu, var, MC, 200);

        EXPECT_FLOAT_EQ(mu[0], 12.3);
        EXPECT_FLOAT_EQ(mu[1], 300);
        EXPECT_FLOAT_EQ(mu[2], 600);

        EXPECT_FLOAT_EQ(var[1], 4500);
        EXPECT_FLOAT_EQ(var[2], 7500);

        EXPECT_FLOAT_EQ(MC[0], 12345);
        EXPECT_FLOAT_EQ(MC[1], 5000);
        EXPECT_FLOAT_EQ(MC[2], 1000);*/
}
