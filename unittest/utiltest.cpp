#include <gtest/gtest.h>
#include "util.h"

TEST(poissonPDF, poissonPDFTest)
{
        EXPECT_DOUBLE_EQ(0.12511003572113349, Util::poissonPDF(10, 10));
        EXPECT_DOUBLE_EQ(4.8646491820674864E-63, Util::poissonPDF(100, 10));
        EXPECT_DOUBLE_EQ(0, Util::poissonPDF(1000, 10));
        EXPECT_DOUBLE_EQ(0, Util::poissonPDF(10000, 10));
}
