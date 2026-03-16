#pragma once
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

class SummaryListener : public Catch::EventListenerBase {
public:
  using EventListenerBase::EventListenerBase;

  void testCaseEnded(Catch::TestCaseStats const& stats) override {
    auto passed = stats.totals.assertions.passed;
    auto failed = stats.totals.assertions.failed;
    auto total = passed + failed;

    if(failed == 0) {
      Catch::cout() << "\033[32m  ✓ " << stats.testInfo->name
                    << ": " << total << " assertion"
                    << (total != 1 ? "s" : "") << " passed\033[0m\n";
    } else {
      Catch::cout() << "\033[31m  ✗ " << stats.testInfo->name
                    << ": " << passed << "/" << total
                    << " passed, " << failed << " FAILED\033[0m\n";
    }
  }
};

CATCH_REGISTER_LISTENER(SummaryListener)
