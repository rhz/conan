#include <conan/utils.hpp>

int main()
{
  using conan::decimal;
  const int num_iter = 1000,
            num_options = 6;
  decimal dist[num_options];
  std::vector<decimal> p(num_options, 1);
  int count = 0;
  double sum = 0,
         percent_sum = 0;
  for (; count < num_options; ++count)
    dist[count] = 0;

  for (count = 0; count < num_iter; ++count)
    ++dist[conan::detail::weighted_choice(p)];

  for (count = 0; count < num_options; ++count)
    sum += dist[count];

  for (count = 0; count < num_options; ++count)
    std::cout << count << '\t' << dist[count] << '\t' << (double)dist[count]/sum << std::endl;

  for (count = 0; count < num_options; ++count)
    percent_sum += (double)dist[count]/sum;

  std::cout << '\t' << sum << '\t' << percent_sum << std::endl;
  return(0);
}
