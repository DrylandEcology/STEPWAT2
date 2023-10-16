#include <gtest/gtest.h>

#include "test_ST_mortality.h"

TEST(ST_Mortality_test, Simulate_prescribed_fire_when_killfreq_startyr_is_0) {
    GrpIndex rg = 0;
    LOG_INFO local_log;
    sw_init_logs(NULL, &LogInfo);

    RGroup = (GroupType **)Mem_Calloc(1, sizeof(GroupType *), nullptr, &local_log);
    RGroup[rg] = (GroupType *)Mem_Calloc(1, sizeof(GroupType), nullptr, &local_log);
    RGroup[rg]->killfreq_startyr = 0;

    simulatePrescribedFire();

    // If killfreq_startyr == 0, do not simulate fire.
    EXPECT_EQ(RGroup[rg]->prescribedfire, 0);

    Mem_Free((void *)RGroup[rg]);
    Mem_Free((void *)RGroup);
}
