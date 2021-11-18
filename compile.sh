echo "Building PochiVM"
cd PochiVM && python3 pochivm-build cmake release && python3 pochivm-build make release && cd ../ && ./update.sh  && \
echo "Building taco" && \
docker run -v `pwd`/:/home/u/taco pochivm-build:latest bash -c 'cd taco/ && cmake . -DBUILD_FLAVOR=RELEASE -GNinja && ninja -j8 && ./bin/taco-tensor_times_vector'
