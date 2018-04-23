#include <gtest/gtest.h>
#include <cstdio>
#include <string>

TEST(Scaffold, CircAPrioriSnap)
{
       /* Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(200, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(1, 3), GaussVal(-100, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(2, 3), GaussVal(-500, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(1), scaffold.getScaffold().size());*/
}

TEST(Scaffold, CircTest)
{
      /*  Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(200, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(1, 3), GaussVal(-100, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(3), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(200, 10), scaffold.getNodeInfo(3).first);*/
}

TEST(Scaffold, ReverseComplPresent)
{
   /*     Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(200, 10));
        scaffold.addObservation(NodePair(1, 4), GaussVal(300, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(2, -3), GaussVal(100, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(2), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(300, 10), scaffold.getNodeInfo(4).first);*/
}

TEST(Scaffold, SnapAPrioriLinkAtRootNode)
{
     /*   Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(200, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(1, 3), GaussVal(1, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(1), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(1).first);*/
}

TEST(Scaffold, SnapPosterioriLinkAtRootNode)
{
     /*   Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(-100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(400, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(1), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(1).first);*/
}

TEST(Scaffold, Linear)
{
    /*    Scaffold scaffold(1);

        EXPECT_EQ(size_t(1), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(1).first);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(-100, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(scaffold.getScaffold().size(), size_t(3));
        EXPECT_EQ(GaussVal(100, 10), scaffold.getNodeInfo(2).first);
        EXPECT_EQ(GaussVal(-100, 10), scaffold.getNodeInfo(3).first);
        EXPECT_EQ(scaffold.getNodeInfo(1).first, GaussVal(0, 0));

        scaffold.addObservation(NodePair(3, 2), GaussVal(200, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(GaussVal(100, 1.0/sqrt(0.015)), scaffold.getNodeInfo(2).first);
        EXPECT_EQ(GaussVal(-100, 1.0/sqrt(0.015)), scaffold.getNodeInfo(3).first);

        scaffold.setRootID(3);
        scaffold.calculateScaffold();

        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(3).first);

        scaffold.setRootID(2);
        scaffold.calculateScaffold();

        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(2).first);*/
}


