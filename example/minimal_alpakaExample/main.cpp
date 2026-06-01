//
// Created by tim on 10.04.26.
//
#include <alpaka/alpaka.hpp>

#include <vector>

static std::vector<std::string> backendList()
{
    auto backends = alpaka::onHost::allBackends(alpaka::onHost::enabledApis, alpaka::exec::enabledExecutors);
    std::vector<std::string> list;
    alpaka::onHost::executeForEachIfHasDevice(
        [&](auto const& backend)
        {
            auto devSelector = alpaka::onHost::makeDeviceSelector(backend[alpaka::object::deviceSpec]);
            auto sampleDevice = devSelector.makeDevice(0);
            std::string backendName;
            backendName += alpaka::onHost::getName(alpaka::getApi(sampleDevice)) + "_";
            backendName += alpaka::onHost::getName(alpaka::getDeviceKind(sampleDevice)) + "_";
            backendName += alpaka::onHost::getName(backend[alpaka::object::exec]);
            list.emplace_back(backendName);
        },
        backends);
    return list;
}

int main(int argc, char const* argv[])
{
    for(auto const& backend : backendList())
    {
        std::cout << backend << std::endl;
    }
}
