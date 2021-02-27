template <typename T>
inline void destroyObject(T &object)
{
  object.~T();
}


template <typename T, typename U>
inline void createObject(T& number, U&& arg)
{
  new (&number) auto(std::forward<U>(arg));
}
