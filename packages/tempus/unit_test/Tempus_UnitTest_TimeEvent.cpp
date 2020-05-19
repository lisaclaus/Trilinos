// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Tempus_TimeEventBase.hpp"
#include "Tempus_TimeEventComposite.hpp"
#include "Tempus_TimeEventRange.hpp"
#include "Tempus_TimeEventRangeIndex.hpp"
#include "Tempus_TimeEventList.hpp"
#include "Tempus_TimeEventListIndex.hpp"

#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <cmath>

static double PI = M_PI;

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventBase, TimeEventBase)
{
  auto te = rcp(new Tempus::TimeEventBase<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventBase");
  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  TEST_COMPARE(te->isTime(0.0), ==, false);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 1.0), ==, false);

  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRange(0, 10), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventRange");
  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_FLOATING_EQUALITY(te->getTimeStart (), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 1);

  TEST_FLOATING_EQUALITY(te->getRelTol(), 1.0e-14, 1.0e-14);
  TEST_COMPARE(te->getLandOnExactly(), ==, true);

  te->setRelTol(0.1);
  TEST_FLOATING_EQUALITY(te->getRelTol(), 0.1, 1.0e-14);
  te->setRelTol(1.0e-14);
  te->setLandOnExactly(false);
  TEST_COMPARE(te->getLandOnExactly(), ==, false);
  te->setLandOnExactly(true);

  // Reset start after stop.
  te->setTimeStart(1.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 1);

  // Reset stop.
  te->setTimeStop(4.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 2);

  // Reset stride.
  te->setTimeStride(0.5);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.5, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 7);

  // Negative stride should be reset to stop_-start_.
  te->setTimeStride(-0.5);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 2);

  // Large stride should be reset to stop_-start_.
  te->setTimeStride(5.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 2);

  // Stride smaller than relative tolerance should be reset to stop_-start_.
  te->setTimeStride(1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 2);

  // Set with time range.
  te->setTimeRange(0.0, PI, 1.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 0.0,  1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 1.0,  1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 4);

  // Test isTime.
  TEST_COMPARE(te->isTime(-10.0e-14), ==, false);   // Around first event.
  TEST_COMPARE(te->isTime( -1.0e-14), ==, true );
  TEST_COMPARE(te->isTime(  0.0    ), ==, true );
  TEST_COMPARE(te->isTime(  1.0e-14), ==, true );
  TEST_COMPARE(te->isTime( 10.0e-14), ==, false);

  TEST_COMPARE(te->isTime(1.0 + -10.0e-14), ==, false);   // Around mid event.
  TEST_COMPARE(te->isTime(1.0 +  -1.0e-14), ==, true );
  TEST_COMPARE(te->isTime(1.0 +   0.0    ), ==, true );
  TEST_COMPARE(te->isTime(1.0 +   1.0e-14), ==, true );
  TEST_COMPARE(te->isTime(1.0 +  10.0e-14), ==, false);

  TEST_COMPARE(te->isTime(3.0 + -10.0e-14), ==, false);   // Around last event.
  TEST_COMPARE(te->isTime(3.0 +  -1.0e-14), ==, true );
  TEST_COMPARE(te->isTime(3.0 +   0.0    ), ==, true );
  TEST_COMPARE(te->isTime(3.0 +   1.0e-14), ==, true );
  TEST_COMPARE(te->isTime(3.0 +  10.0e-14), ==, false);

  // Test timeToNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-10.0e-14),  1.0e-13, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( -1.0e-14),  1.0e-14, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(  0.0    ),  0.0    , 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(  1.0e-14), -1.0e-14, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( 10.0e-14), 1.0-1.0e-13, 1.0e-14);

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+ -10.0e-14),  1.0e-13, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+  -1.0e-14),  1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+   0.0    ),  0.0    , 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+   1.0e-14), -1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+  10.0e-14), 1.0-1.0e-13, 1.0e-14);

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+ -10.0e-14),  1.0e-13, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+  -1.0e-14),  1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+   0.0    ),  0.0    , 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+   1.0e-14), -1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+  10.0e-14), -1.0e-13, 1.0e-02);

  // Test timeOfNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-10.0e-14), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( -1.0e-14), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(  0.0    ), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(  1.0e-14), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( 10.0e-14), 1.0, 1.0e-14);

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+ -10.0e-14), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+  -1.0e-14), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+   0.0    ), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+   1.0e-14), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+  10.0e-14), 2.0, 1.0e-14);

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+ -10.0e-14), 3.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+  -1.0e-14), 3.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+   0.0    ), 3.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+   1.0e-14), 3.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+  10.0e-14), 3.0, 1.0e-14);

  // Test eventInRange.
  //   Right end.
  TEST_COMPARE(te->eventInRange(-1.0, -10.0e-14), ==, false);   // Around first event.
  TEST_COMPARE(te->eventInRange(-1.0,  -1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(-1.0,   0.0    ), ==, true );
  TEST_COMPARE(te->eventInRange(-1.0,   1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(-1.0,  10.0e-14), ==, true );

  TEST_COMPARE(te->eventInRange(0.5, 1.0 + -10.0e-14), ==, false);   // Around mid event.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 +  -1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(0.5, 1.0 +   0.0    ), ==, true );
  TEST_COMPARE(te->eventInRange(0.5, 1.0 +   1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(0.5, 1.0 +  10.0e-14), ==, true );

  TEST_COMPARE(te->eventInRange(2.5, 3.0 + -10.0e-14), ==, false);   // Around last event.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 +  -1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(2.5, 3.0 +   0.0    ), ==, true );
  TEST_COMPARE(te->eventInRange(2.5, 3.0 +   1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(2.5, 3.0 +  10.0e-14), ==, true );

  //   Left end.
  TEST_COMPARE(te->eventInRange(-10.0e-14, 0.5), ==, true );   // Around first event.
  TEST_COMPARE(te->eventInRange( -1.0e-14, 0.5), ==, true );
  TEST_COMPARE(te->eventInRange(  0.0    , 0.5), ==, true );
  TEST_COMPARE(te->eventInRange(  1.0e-14, 0.5), ==, true );
  TEST_COMPARE(te->eventInRange( 10.0e-14, 0.5), ==, false );

  TEST_COMPARE(te->eventInRange(1.0 + -10.0e-14, 1.5), ==, true );   // Around mid event.
  TEST_COMPARE(te->eventInRange(1.0 +  -1.0e-14, 1.5), ==, true );
  TEST_COMPARE(te->eventInRange(1.0 +   0.0    , 1.5), ==, true );
  TEST_COMPARE(te->eventInRange(1.0 +   1.0e-14, 1.5), ==, true );
  TEST_COMPARE(te->eventInRange(1.0 +  10.0e-14, 1.5), ==, false);

  TEST_COMPARE(te->eventInRange(3.0 + -10.0e-14, 4.0), ==, true );   // Around last event.
  TEST_COMPARE(te->eventInRange(3.0 +  -1.0e-14, 4.0), ==, true );
  TEST_COMPARE(te->eventInRange(3.0 +   0.0    , 4.0), ==, true );
  TEST_COMPARE(te->eventInRange(3.0 +   1.0e-14, 4.0), ==, true );
  TEST_COMPARE(te->eventInRange(3.0 +  10.0e-14, 4.0), ==, false);

  // Check base class defaults (functions not implemented in TimeEventRange).
  TEST_COMPARE(te->isIndex(1), ==, false);
  TEST_COMPARE(te->indexToNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRangeIndex(1,4), ==, false);

}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Construction_Stride)
{
  auto te = rcp(new Tempus::TimeEventRange<double>(
                "TestName", 0.0, PI, 1.0, 1.0e-14, true));

  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 0.0 , 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 1.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 4);

  auto teRange2 = rcp(new Tempus::TimeEventRange<double>(
      "teRange2", -PI/2.0, PI/2.0, PI/4.0, 1.0e-14, true));

  TEST_FLOATING_EQUALITY(teRange2->timeToNextEvent(0.1), PI/4.0-0.1, 1.0e-14);

}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Construction_NumEvents)
{
  auto te = rcp(new Tempus::TimeEventRange<double>(
                "TestName", 0.0, PI, 5, 1.0e-14, true));

  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 0.0,      1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), PI,     1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), PI/4.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 5);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventList<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventList");
  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  TEST_COMPARE(te->getTimeList().size(), ==, 0);
  TEST_FLOATING_EQUALITY(te->getRelTol(), 1.0e-14, 1.0e-14);
  TEST_COMPARE(te->getLandOnExactly(), ==, true);

  te->setRelTol(0.1);
  TEST_FLOATING_EQUALITY(te->getRelTol(), 0.1, 1.0e-14);
  te->setRelTol(1.0e-14);
  te->setLandOnExactly(false);
  TEST_COMPARE(te->getLandOnExactly(), ==, false);
  te->setLandOnExactly(true);


  // Test isTime with zero elements.
  TEST_COMPARE(te->isTime(te->getDefaultTime()*(1.0 + -10.0e-14)), ==, false);
  TEST_COMPARE(te->isTime(te->getDefaultTime()*(1.0 +  -1.0e-14)), ==, true );
  TEST_COMPARE(te->isTime(te->getDefaultTime()*(1.0 +   0.0    )), ==, true );
  TEST_COMPARE(te->isTime(te->getDefaultTime()*(1.0 +   1.0e-14)), ==, true );
  TEST_COMPARE(te->isTime(te->getDefaultTime()*(1.0 +  10.0e-14)), ==, false);

  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_COMPARE(te->eventInRange(2.0, -1.0), ==, false);

  // Test isTime with one element.
  te->addTime(0.0);
  TEST_COMPARE(te->isTime(-10.0e-14), ==, false);
  TEST_COMPARE(te->isTime( -1.0e-14), ==, true );
  TEST_COMPARE(te->isTime(  0.0    ), ==, true );
  TEST_COMPARE(te->isTime(  1.0e-14), ==, true );
  TEST_COMPARE(te->isTime( 10.0e-14), ==, false);

  // Test isTime with two elements.
  te->addTime(PI);
  TEST_COMPARE(te->isTime(PI + -10.0e-14), ==, false);
  TEST_COMPARE(te->isTime(PI +  -1.0e-14), ==, true );
  TEST_COMPARE(te->isTime(PI +   0.0    ), ==, true );
  TEST_COMPARE(te->isTime(PI +   1.0e-14), ==, true );
  TEST_COMPARE(te->isTime(PI +  10.0e-14), ==, false);

  // Test addTime.
  te->addTime(-1.0);
  te->addTime( 2.0);
  te->addTime( 5.0);
  auto testList = te->getTimeList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_FLOATING_EQUALITY(testList[0], -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[1],  0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[2],  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[3], PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[4],  5.0, 1.0e-14);

  // Test timeToNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 + -10.0e-14),  1.0e-13, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 +  -1.0e-14),  1.0e-14, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 +   0.0    ),  0.0    , 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 +   1.0e-14), -1.0e-14, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 +  10.0e-14), 1.0-1.0e-13, 1.0e-14);

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI + -10.0e-14),  1.0e-13, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI +  -1.0e-14),  1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI +   0.0    ),  0.0    , 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI +   1.0e-14), -1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI +  10.0e-14), 5.0-PI-1.0e-13, 1.0e-14);

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+ -10.0e-14),  1.0e-13, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+  -1.0e-14),  1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+   0.0    ),  0.0    , 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+   1.0e-14), -1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+  10.0e-14), -1.0e-13, 1.0e-02);

  // Test timeOfNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 + -10.0e-14), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 +  -1.0e-14), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 +   0.0    ), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 +   1.0e-14), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 +  10.0e-14),  0.0, 1.0e-14);

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+ -10.0e-14),  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+  -1.0e-14),  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+   0.0    ),  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+   1.0e-14),  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+  10.0e-14), PI, 1.0e-14);

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+ -10.0e-14), 5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+  -1.0e-14), 5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+   0.0    ), 5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+   1.0e-14), 5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+  10.0e-14), 5.0, 1.0e-14);

  // Test eventInRange.
  //   Right end.
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 + -10.0e-14), ==, false);   // Around first event.
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 +  -1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 +   0.0    ), ==, true );
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 +   1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 +  10.0e-14), ==, true );

  TEST_COMPARE(te->eventInRange(3.0, PI + -10.0e-14), ==, false);   // Around mid event.
  TEST_COMPARE(te->eventInRange(3.0, PI +  -1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(3.0, PI +   0.0    ), ==, true );
  TEST_COMPARE(te->eventInRange(3.0, PI +   1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(3.0, PI +  10.0e-14), ==, true );

  TEST_COMPARE(te->eventInRange(4.5, 5.0 + -10.0e-14), ==, false);   // Around last event.
  TEST_COMPARE(te->eventInRange(4.5, 5.0 +  -1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(4.5, 5.0 +   0.0    ), ==, true );
  TEST_COMPARE(te->eventInRange(4.5, 5.0 +   1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(4.5, 5.0 +  10.0e-14), ==, true );

  //   Left end.
  TEST_COMPARE(te->eventInRange(-1.0 + -10.0e-14, -0.5), ==, true );   // Around first event.
  TEST_COMPARE(te->eventInRange(-1.0 +  -1.0e-14, -0.5), ==, true );
  TEST_COMPARE(te->eventInRange(-1.0 +   0.0    , -0.5), ==, true );
  TEST_COMPARE(te->eventInRange(-1.0 +   1.0e-14, -0.5), ==, true );
  TEST_COMPARE(te->eventInRange(-1.0 +  10.0e-14, -0.5), ==, false );

  TEST_COMPARE(te->eventInRange(PI + -10.0e-14, 3.5), ==, true );   // Around mid event.
  TEST_COMPARE(te->eventInRange(PI +  -1.0e-14, 3.5), ==, true );
  TEST_COMPARE(te->eventInRange(PI +   0.0    , 3.5), ==, true );
  TEST_COMPARE(te->eventInRange(PI +   1.0e-14, 3.5), ==, true );
  TEST_COMPARE(te->eventInRange(PI +  10.0e-14, 3.5), ==, false);

  TEST_COMPARE(te->eventInRange(5.0 + -10.0e-14, 6.0), ==, true );   // Around last event.
  TEST_COMPARE(te->eventInRange(5.0 +  -1.0e-14, 6.0), ==, true );
  TEST_COMPARE(te->eventInRange(5.0 +   0.0    , 6.0), ==, true );
  TEST_COMPARE(te->eventInRange(5.0 +   1.0e-14, 6.0), ==, true );
  TEST_COMPARE(te->eventInRange(5.0 +  10.0e-14, 6.0), ==, false);

  // Check base class defaults (functions not implemented in TimeEventList).
  TEST_COMPARE(te->isIndex(1), ==, false);
  TEST_COMPARE(te->indexToNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRangeIndex(1,4), ==, false);

}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, Construction)
{
  std::vector<double> testVector;
  testVector.push_back(-1.0);
  testVector.push_back( 0.0);
  testVector.push_back( 5.0);
  testVector.push_back( 2.0);
  testVector.push_back(PI);

  auto te = rcp(new Tempus::TimeEventList<double>(
                "TestName", testVector, 1.0e-14, true));

  TEST_COMPARE(te->getName(), ==, "TestName");
  TEST_FLOATING_EQUALITY(te->getRelTol(), 1.0e-14, 1.0e-14);
  TEST_COMPARE(te->getLandOnExactly(), ==, true);

  auto testList = te->getTimeList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_FLOATING_EQUALITY(testList[0], -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[1],  0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[2],  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[3], PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[4],  5.0, 1.0e-14);

  // Test that two events within relative tolerance are added or not.
  te->addTime( 2.0 + 1.0e-14);
  TEST_COMPARE(te->getTimeList().size(), ==, 5);
  te->addTime( 2.0 + 1.0e-13);
  TEST_COMPARE(te->getTimeList().size(), ==, 6);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventRangeIndex");
  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_COMPARE(te->getIndexStart (), ==, 0);
  TEST_COMPARE(te->getIndexStop  (), ==, 0);
  TEST_COMPARE(te->getIndexStride(), ==, 1);
  TEST_COMPARE(te->getNumEvents  (), ==, 1);

  // Reset start after stop.
  te->setIndexStart(1);
  TEST_COMPARE(te->getIndexStart (), ==, 1);
  TEST_COMPARE(te->getIndexStop  (), ==, 1);
  TEST_COMPARE(te->getIndexStride(), ==, 1);
  TEST_COMPARE(te->getNumEvents  (), ==, 1);

  // Reset stop.
  te->setIndexStop(5);
  TEST_COMPARE(te->getIndexStart (), ==, 1);
  TEST_COMPARE(te->getIndexStop  (), ==, 5);
  TEST_COMPARE(te->getIndexStride(), ==, 1);
  TEST_COMPARE(te->getNumEvents  (), ==, 5);

  // Reset stride.
  te->setIndexStride(2);
  TEST_COMPARE(te->getIndexStart (), ==, 1);
  TEST_COMPARE(te->getIndexStop  (), ==, 5);
  TEST_COMPARE(te->getIndexStride(), ==, 2);
  TEST_COMPARE(te->getNumEvents  (), ==, 3);

  // Negative stride should be reset to stop_-start_.
  te->setIndexStride(-1);
  TEST_COMPARE(te->getIndexStart (), ==, 1);
  TEST_COMPARE(te->getIndexStop  (), ==, 5);
  TEST_COMPARE(te->getIndexStride(), ==, 1);
  TEST_COMPARE(te->getNumEvents  (), ==, 5);

  // Large stride should be reset to stop_-start_.
  te->setIndexStride(5);
  TEST_COMPARE(te->getIndexStart (), ==, 1);
  TEST_COMPARE(te->getIndexStop  (), ==, 5);
  TEST_COMPARE(te->getIndexStride(), ==, 4);
  TEST_COMPARE(te->getNumEvents  (), ==, 2);

  // Set with index range.
  te->setIndexRange(-5, 5, 3);
  TEST_COMPARE(te->getIndexStart (), ==, -5);
  TEST_COMPARE(te->getIndexStop  (), ==,  5);
  TEST_COMPARE(te->getIndexStride(), ==,  3);
  TEST_COMPARE(te->getNumEvents  (), ==,  4);

  // Test isIndex.
  TEST_COMPARE(te->isIndex(-6), ==, false);   // Around first event.
  TEST_COMPARE(te->isIndex(-5), ==, true );
  TEST_COMPARE(te->isIndex(-4), ==, false);

  TEST_COMPARE(te->isIndex(0), ==, false);    // Around mid event.
  TEST_COMPARE(te->isIndex(1), ==, true );
  TEST_COMPARE(te->isIndex(2), ==, false);

  TEST_COMPARE(te->isIndex(3), ==, false);    // Around last event.
  TEST_COMPARE(te->isIndex(4), ==, true );
  TEST_COMPARE(te->isIndex(5), ==, false);

  // Test indexToNextEvent.
  TEST_COMPARE(te->indexToNextEvent(-9), ==, 4); //   Around first event.
  TEST_COMPARE(te->indexToNextEvent(-5), ==, 0);
  TEST_COMPARE(te->indexToNextEvent(-4), ==, 2);

  TEST_COMPARE(te->indexToNextEvent(-1), ==, 2); //   Around mid event.
  TEST_COMPARE(te->indexToNextEvent( 1), ==, 0);
  TEST_COMPARE(te->indexToNextEvent( 3), ==, 1);

  TEST_COMPARE(te->indexToNextEvent( 2), ==, 2); //   Around last event.
  TEST_COMPARE(te->indexToNextEvent( 4), ==, 0);
  TEST_COMPARE(te->indexToNextEvent( 8), ==,-4);

  // Test indexOfNextEvent.
  TEST_COMPARE(te->indexOfNextEvent(-9), ==, -5); //   Around first event.
  TEST_COMPARE(te->indexOfNextEvent(-5), ==, -5);
  TEST_COMPARE(te->indexOfNextEvent(-4), ==, -2);

  TEST_COMPARE(te->indexOfNextEvent(-1), ==, 1); //   Around mid event.
  TEST_COMPARE(te->indexOfNextEvent( 1), ==, 1);
  TEST_COMPARE(te->indexOfNextEvent( 3), ==, 4);

  TEST_COMPARE(te->indexOfNextEvent( 2), ==, 4); //   Around last event.
  TEST_COMPARE(te->indexOfNextEvent( 4), ==, 4);
  TEST_COMPARE(te->indexOfNextEvent( 8), ==, 4);

  // Test eventInRangeIndex.
  //   Right end.
  TEST_COMPARE(te->eventInRangeIndex(-9, -6), ==, false);   // Around first event.
  TEST_COMPARE(te->eventInRangeIndex(-9, -5), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-9, -4), ==, true );

  TEST_COMPARE(te->eventInRangeIndex(-1, 0), ==, false);   // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-1, 1), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-1, 2), ==, true );

  TEST_COMPARE(te->eventInRangeIndex(2, 3), ==, false);   // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(2, 4), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(2, 5), ==, true );

  //   Left end.
  TEST_COMPARE(te->eventInRangeIndex(-6.0, -3), ==, true );   // Around first event.
  TEST_COMPARE(te->eventInRangeIndex(-5.0, -3), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-4.0, -3), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(-3, 0), ==, true );   // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-2, 0), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-1, 0), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(3, 8), ==, true );   // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(4, 8), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(5, 8), ==, false);

  // Check base class defaults (functions not implemented in TimeEventRangeIndex).
  TEST_COMPARE(te->isTime(1.0), ==, false);
  TEST_COMPARE(te->timeToNextEvent(1.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->timeOfNextEvent(1.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->eventInRange(1.0, 4.0), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, Construction)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>(
                "TestName", -1, 10, 2));

  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_COMPARE(te->getIndexStart (), ==, -1);
  TEST_COMPARE(te->getIndexStop  (), ==, 10);
  TEST_COMPARE(te->getIndexStride(), ==,  2);
  TEST_COMPARE(te->getNumEvents  (), ==,  6);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventListIndex<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventListIndex");
  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  TEST_COMPARE(te->getIndexList().size(), ==, 0);


  // Test isTime with zero elements.
  TEST_COMPARE(te->isIndex(-1), ==, false);
  TEST_COMPARE(te->isIndex(te->getDefaultIndex()), ==, true );
  TEST_COMPARE(te->isIndex( 1), ==, false);

  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRangeIndex(2, -1), ==, false);

  // Test isTime with one element.
  te->addIndex(-5);
  TEST_COMPARE(te->isIndex(-6), ==, false);
  TEST_COMPARE(te->isIndex(-5), ==, true );
  TEST_COMPARE(te->isIndex(-4), ==, false);

  // Test isTime with two elements.
  te->addIndex(1);
  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->isIndex(1), ==, true );
  TEST_COMPARE(te->isIndex(2), ==, false);

  // Test addIndex.
  te->addIndex(-2);
  te->addIndex( 4);
  te->addIndex(-9);
  auto testList = te->getIndexList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_COMPARE(testList[0], ==, -9);
  TEST_COMPARE(testList[1], ==, -5);
  TEST_COMPARE(testList[2], ==, -2);
  TEST_COMPARE(testList[3], ==,  1);
  TEST_COMPARE(testList[4], ==,  4);

  // Test indexToNextEvent.
  //   Around first event.
  TEST_COMPARE(te->indexToNextEvent(-12), ==,  3);
  TEST_COMPARE(te->indexToNextEvent( -9), ==,  0);
  TEST_COMPARE(te->indexToNextEvent( -8), ==,  3);

  //   Around mid event.
  TEST_COMPARE(te->indexToNextEvent(-4), ==,  2);
  TEST_COMPARE(te->indexToNextEvent(-2), ==,  0);
  TEST_COMPARE(te->indexToNextEvent( 0), ==,  1);

  //   Around last event.
  TEST_COMPARE(te->indexToNextEvent(2), ==,  2);
  TEST_COMPARE(te->indexToNextEvent(4), ==,  0);
  TEST_COMPARE(te->indexToNextEvent(9), ==, -5);

  // Test indexOfNextEvent.
  //   Around first event.
  TEST_COMPARE(te->indexOfNextEvent(-12), ==, -9);
  TEST_COMPARE(te->indexOfNextEvent( -9), ==, -9);
  TEST_COMPARE(te->indexOfNextEvent( -8), ==, -5);

  //   Around mid event.
  TEST_COMPARE(te->indexOfNextEvent(-4), ==, -2);
  TEST_COMPARE(te->indexOfNextEvent(-2), ==, -2);
  TEST_COMPARE(te->indexOfNextEvent( 0), ==,  1);

  //   Around last event.
  TEST_COMPARE(te->indexOfNextEvent(2), ==,  4);
  TEST_COMPARE(te->indexOfNextEvent(4), ==,  4);
  TEST_COMPARE(te->indexOfNextEvent(9), ==,  4);

  // Test eventInRangeIndex.
  //   Right end.
  TEST_COMPARE(te->eventInRangeIndex(-12.0, -10), ==, false);   // Around first event.
  TEST_COMPARE(te->eventInRangeIndex(-12.0,  -9), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-12.0,  -8), ==, true );

  TEST_COMPARE(te->eventInRangeIndex(-4, -3), ==, false);   // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-4, -2), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-4, -1), ==, true );

  TEST_COMPARE(te->eventInRangeIndex(3, 3), ==, false);   // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(3, 4), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(3, 6), ==, true );

  //   Left end.
  TEST_COMPARE(te->eventInRangeIndex(-12, -7), ==, true );   // Around first event.
  TEST_COMPARE(te->eventInRangeIndex( -9, -7), ==, true );
  TEST_COMPARE(te->eventInRangeIndex( -8, -7), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(-3, 0), ==, true );   // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-2, 0), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-1, 0), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(3, 8), ==, true );   // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(4, 8), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(5, 8), ==, false);

  // Check base class defaults (functions not implemented in TimeEventListIndex).
  TEST_COMPARE(te->isTime(1.0), ==, false);
  TEST_COMPARE(te->timeToNextEvent(1.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->timeOfNextEvent(1.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->eventInRange(1.0, 4.0), ==, false);

}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, Construction)
{
  std::vector<int> testVector;
  testVector.push_back(-2);
  testVector.push_back( 0);
  testVector.push_back( 7);
  testVector.push_back( 3);
  testVector.push_back(-5);

  auto te = rcp(new Tempus::TimeEventListIndex<double>(
                "TestName", testVector));

  TEST_COMPARE(te->getName(), ==, "TestName");

  auto testList = te->getIndexList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_COMPARE(testList[0], ==, -5);
  TEST_COMPARE(testList[1], ==, -2);
  TEST_COMPARE(testList[2], ==,  0);
  TEST_COMPARE(testList[3], ==,  3);
  TEST_COMPARE(testList[4], ==,  7);

  // Test adding a duplicate event index.
  te->addIndex(3);
  TEST_COMPARE(te->getIndexList().size(), ==, 5);
  te->addIndex(1);
  TEST_COMPARE(te->getIndexList().size(), ==, 6);
}



// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, Construction)
{
  // TimeEventRanges for testing.
  auto teRange1 = rcp(new Tempus::TimeEventRange<double>(
    "teRange1", 0.0, PI, 1.0, 1.0e-14, true));

  auto teRange2 = rcp(new Tempus::TimeEventRange<double>(
    "teRange2", -PI/2.0, PI/2.0, PI/4.0, 1.0e-14, true));

  auto teRange3 = rcp(new Tempus::TimeEventRange<double>(
    "teRange3", 4.0, 10.0, 4.0, 1.0e-14, true));

  // TimeEventLists for testing.
  std::vector<double> testList1;
  testList1.push_back(-1.0);
  testList1.push_back( 0.0);
  testList1.push_back( 5.0);
  testList1.push_back( 2.0);
  testList1.push_back(PI);

  auto teList1 = rcp(new Tempus::TimeEventList<double>(
                     "teList1", testList1, 1.0e-14, true));

  std::vector<double> testList2;
  testList2.push_back(-0.5);
  testList2.push_back( 1.25);
  testList2.push_back( 4.95);
  testList2.push_back(12.34);

  auto teList2 = rcp(new Tempus::TimeEventList<double>(
                     "teList2", testList2, 1.0e-14, true));

  std::vector<double> testList3;
  testList3.push_back(-5.0);
  testList3.push_back(-PI);

  auto teList3 = rcp(new Tempus::TimeEventList<double>(
                     "teList3", testList3, 1.0e-14, true));


  auto teRangeIndex1 = rcp(new Tempus::TimeEventRangeIndex<double>(
                           "teRangeIndex1", -1, 10, 3));

  auto teRangeIndex2 = rcp(new Tempus::TimeEventRangeIndex<double>(
                           "teRangeIndex2", -5, 8, 4));

  auto teRangeIndex3 = rcp(new Tempus::TimeEventRangeIndex<double>(
                           "teRangeIndex2", 10, 17, 5));

  std::vector<int> testListIndex1;
  testListIndex1.push_back(-2);
  testListIndex1.push_back( 0);
  testListIndex1.push_back( 7);
  testListIndex1.push_back( 3);
  testListIndex1.push_back(-5);

  auto teListIndex1 = rcp(new Tempus::TimeEventListIndex<double>(
                          "teListIndex1", testListIndex1));

  std::vector<int> testListIndex2;
  testListIndex2.push_back( 2);

  auto teListIndex2 = rcp(new Tempus::TimeEventListIndex<double>(
                          "teListIndex2", testListIndex2));

  std::vector<int> testListIndex3;
  testListIndex3.push_back(14);
  testListIndex3.push_back( 9);

  auto teListIndex3 = rcp(new Tempus::TimeEventListIndex<double>(
                          "teListIndex3", testListIndex3));

  // ---------------------------------------------------------------------------

  // Check with zero TimeEvents.
  auto te = rcp(new Tempus::TimeEventComposite<double>());

  TEST_COMPARE(te->isTime(0.0), ==, false);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 1.0), ==, false);

  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRange(0, 10), ==, false);

  // Check with one TimeEventRange.
  te->clearTimeEvents();
  te->addTimeEvent(teRange1);

  TEST_COMPARE(te->isTime(0.0), ==, true );
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.0), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.1), 2.0, 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 1.0), ==, true );

  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRange(0, 10), ==, true );

  // Check with two overlapping TimeEventRanges.
  te->clearTimeEvents();
  te->addTimeEvent(teRange1);
  te->addTimeEvent(teRange2);
  //te->describe();

  TEST_COMPARE(te->isTime(0.0), ==, true );
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.1), PI/4.0-0.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.1), PI/2.0, 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 1.0), ==, true );

  // Check with two separated TimeEventRanges.
  te->clearTimeEvents();
  te->addTimeEvent(teRange3);
  te->addTimeEvent(teRange2);
  //te->describe();

  TEST_COMPARE(te->isTime(4.0), ==, true );
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.1), PI/2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0),    4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0),    8.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(10.0),   8.0, 1.0e-14);
  TEST_COMPARE(te->eventInRange(-3.0, -2.0), ==, false);
  TEST_COMPARE(te->eventInRange(-1.0,  0.0), ==, true );
  TEST_COMPARE(te->eventInRange( 1.0,  2.0), ==, true );
  TEST_COMPARE(te->eventInRange( 2.0,  3.0), ==, false);
  TEST_COMPARE(te->eventInRange( 5.0,  7.0), ==, false);
  TEST_COMPARE(te->eventInRange( 7.0,  9.0), ==, true );
  TEST_COMPARE(te->eventInRange( 9.0, 11.0), ==, false);

  // ---------------------------------------------------------------------------

  // Check with one TimeEventList.
  te->clearTimeEvents();
  te->addTimeEvent(teList1);

  TEST_COMPARE(te->isTime(2.0), ==, true );
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.5), 1.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-2.0), -1.0, 1.0e-14);
  TEST_COMPARE(te->eventInRange(4.99, 10.0), ==, true );

  // Check with two overlapping TimeEventLists.
  te->clearTimeEvents();
  te->addTimeEvent(teList1);
  te->addTimeEvent(teList2);
  //te->describe();

  TEST_COMPARE(te->isTime(1.25), ==, true );
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(4.0), 0.95, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(6.5), 12.34, 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.1, 1.0), ==, false);

  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRange(0, 10), ==, true );

  // Check with two separated TimeEventLists.
  te->clearTimeEvents();
  te->addTimeEvent(teList3);
  te->addTimeEvent(teList2);
  //te->describe();

  TEST_COMPARE(te->isTime(4.0), ==, false);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-8.9),  -5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-0.3),  1.25, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-4.0),   -PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(20.0), 12.34, 1.0e-14);
  TEST_COMPARE(te->eventInRange(-6.0, -4.0), ==, true );
  TEST_COMPARE(te->eventInRange(-3.0,  0.0), ==, true );
  TEST_COMPARE(te->eventInRange( 2.0,  3.0), ==, false);
  TEST_COMPARE(te->eventInRange( 4.9,  5.1), ==, true );
  TEST_COMPARE(te->eventInRange(12.0, 12.4), ==, true );
  TEST_COMPARE(te->eventInRange(14.0, 15.0), ==, false);

  // ---------------------------------------------------------------------------

  // Check with one TimeEventRangeIndex.
  te->clearTimeEvents();
  te->addTimeEvent(teRangeIndex1);

  TEST_COMPARE(te->isIndex(5), ==, true );
  TEST_COMPARE(te->indexToNextEvent(3), ==, 2);
  TEST_COMPARE(te->indexOfNextEvent(3), ==, 5);
  TEST_COMPARE(te->eventInRangeIndex(3, 9), ==, true );

  TEST_COMPARE(te->isTime(0.0), ==, false);
  TEST_COMPARE(te->timeToNextEvent(0.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->timeOfNextEvent(0.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->eventInRange(0.0, 10.0), ==, false);

  // Check with two overlapping TimeEventRangeIndices.
  te->clearTimeEvents();
  te->addTimeEvent(teRangeIndex1);
  te->addTimeEvent(teRangeIndex2);
  //te->describe();

  TEST_COMPARE(te->isIndex(-1), ==, true );
  TEST_COMPARE(te->indexToNextEvent(-2), ==, 1);
  TEST_COMPARE(te->indexOfNextEvent( 2), ==, 2);
  TEST_COMPARE(te->eventInRangeIndex(0, 1), ==, false);

  // Check with two separated TimeEventRangesIndices.
  te->clearTimeEvents();
  te->addTimeEvent(teRangeIndex3);
  te->addTimeEvent(teRangeIndex2);
  //te->describe();

  TEST_COMPARE(te->isIndex(15), ==, true );
  TEST_COMPARE(te->indexOfNextEvent( 9), ==, 10);
  TEST_COMPARE(te->indexOfNextEvent(-6), ==, -5);
  TEST_COMPARE(te->indexOfNextEvent( 6), ==,  7);
  TEST_COMPARE(te->indexOfNextEvent(16), ==, 15);
  TEST_COMPARE(te->eventInRangeIndex(-3, -2), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(-1,  0), ==, true );
  TEST_COMPARE(te->eventInRangeIndex( 1,  2), ==, false);
  TEST_COMPARE(te->eventInRangeIndex( 7,  8), ==, true );
  TEST_COMPARE(te->eventInRangeIndex( 8,  9), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(10, 13), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(14, 20), ==, true );

  // ---------------------------------------------------------------------------

  // Check with one TimeEventListIndex.
  te->clearTimeEvents();
  te->addTimeEvent(teListIndex1);

  TEST_COMPARE(te->isIndex(3), ==, true );
  TEST_COMPARE(te->indexToNextEvent(1), ==, 2);
  TEST_COMPARE(te->indexOfNextEvent(4), ==, 7);
  TEST_COMPARE(te->eventInRangeIndex(1, 3), ==, true );

  // Check with two overlapping TimeEventListIndices.
  te->clearTimeEvents();
  te->addTimeEvent(teListIndex1);
  te->addTimeEvent(teListIndex2);
  //te->describe();

  TEST_COMPARE(te->isIndex(2), ==, true );
  TEST_COMPARE(te->indexToNextEvent(0), ==, 0);
  TEST_COMPARE(te->indexOfNextEvent(1), ==, 2);
  TEST_COMPARE(te->eventInRangeIndex(-1, 3), ==, true );

  TEST_COMPARE(te->isTime(0.0), ==, false);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(10.0), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-3.0), te->getDefaultTime(), 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 10.0), ==, false);

  // Check with two separated TimeEventListIndices.
  te->clearTimeEvents();
  te->addTimeEvent(teListIndex3);
  te->addTimeEvent(teListIndex2);
  //te->describe();

  TEST_COMPARE(te->isIndex(14), ==, true );
  TEST_COMPARE(te->indexOfNextEvent(2), ==, 2);
  TEST_COMPARE(te->indexOfNextEvent(5), ==, 9);
  TEST_COMPARE(te->indexOfNextEvent(19), ==, 14);
  TEST_COMPARE(te->eventInRangeIndex( 0,  1), ==, false);
  TEST_COMPARE(te->eventInRangeIndex( 3, 10), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(15, 20), ==, false);

  // ---------------------------------------------------------------------------

  // Check with one of everything.
  te->clearTimeEvents();
  te->addTimeEvent(teRange1);
  te->addTimeEvent(teList1);
  te->addTimeEvent(teRangeIndex1);
  te->addTimeEvent(teListIndex1);

  TEST_COMPARE(te->isTime (3.0), ==, true );
  TEST_COMPARE(te->isTime (2.0), ==, true );
  TEST_COMPARE(te->isIndex(  2), ==, true );
  TEST_COMPARE(te->isIndex(  3), ==, true );

  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-2.5),  1.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( 0.5),  0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( 4.5),  0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( 7.5), -2.5, 1.0e-14);

  TEST_COMPARE(te->indexToNextEvent(-6), ==,  1);
  TEST_COMPARE(te->indexToNextEvent( 1), ==,  1);
  TEST_COMPARE(te->indexToNextEvent( 7), ==,  0);
  TEST_COMPARE(te->indexToNextEvent( 9), ==, -1);

  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( -PI), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-0.5),  0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( 2.5),  3.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( 7.5),  5.0, 1.0e-14);

  TEST_COMPARE(te->indexOfNextEvent(-6), ==, -5);
  TEST_COMPARE(te->indexOfNextEvent( 1), ==,  2);
  TEST_COMPARE(te->indexOfNextEvent( 7), ==,  7);
  TEST_COMPARE(te->indexOfNextEvent( 9), ==,  8);

  TEST_COMPARE(te->eventInRange(-5.0, -2.0), ==, false);
  TEST_COMPARE(te->eventInRange(-2.0, -0.5), ==, true );
  TEST_COMPARE(te->eventInRange( 1.2,  1.8), ==, false);
  TEST_COMPARE(te->eventInRange( 3.1,  4.0), ==, true );
  TEST_COMPARE(te->eventInRange( 4.5,  6.0), ==, true );

  TEST_COMPARE(te->eventInRangeIndex(-8, -6), ==, false);
  TEST_COMPARE(te->eventInRangeIndex( 1,  1), ==, false);
  TEST_COMPARE(te->eventInRangeIndex( 5,  7), ==, true );
  TEST_COMPARE(te->eventInRangeIndex( 8, 10), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(12, 14), ==, false);

}


} // namespace Tempus_Test
