int fib(const int&n, std::ostream* pstream__) {
  if (n <= 0) {
    stringstream errmsg;
    errmsg << "n must be positive";
    throw std::domain_error(errmsg.str());
  }
  return n <= 1 ? 1 : fib(n - 1, 0) + fib(n - 2, 0);
}
