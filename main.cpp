#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include <sstream>

using std::vector;
using std::cerr;
using std::ostream;
using Number = float;


template <typename T>
static ostream& operator<<(ostream& s, const vector<T> &v)
{
  s << "[";

  for (const auto &x : v) {
    s << " " << x;
  }

  s << " ]";
  return s;
}


namespace {
struct Array;
}


static bool operator==(const Array &a, const Array &b);


template <typename T>
static void destroyObject(T &object)
{
  object.~T();
}


namespace {
struct Value {
  enum class Type {
    number, character, array_ptr
  };

  friend ostream& operator<<(ostream &stream, Type type)
  {
    switch (type) {
      case Type::number:
        stream << "number";
        break;
      case Type::character:
        assert(false);
      case Type::array_ptr:
        stream << "array_ptr";
        break;
    }

    return stream;
  }

  Type type;

  union {
    Number number;
    char character;
    std::unique_ptr<Array> array_ptr;
  };

  Value()
  : type(Type::array_ptr), array_ptr(nullptr)
  {
  }

  Value(Number arg)
  : type(Type::number), number(arg)
  {
  }

  Value(int arg)
  : type(Type::number), number(arg)
  {
  }

  Value(char arg)
  : type(Type::character), character(arg)
  {
  }

  Value(Array arg);

  Value(Value &&arg)
  : type(arg.type)
  {
    createFrom(std::move(arg));
  }

  explicit Value(const Value &arg)
  : type(arg.type)
  {
    createFrom(arg);
  }

  bool isArray() const { return type == Type::array_ptr; }
  bool isNumber() const { return type == Type::number; }
  bool isCharacter() const { return type == Type::character; }

  Number asNumber() const
  {
    assert(type == Type::number);
    return number;
  }

  char asCharacter() const
  {
    assert(type == Type::character);
    return character;
  }

  const Array &asArray() const
  {
    assert(type == Type::array_ptr);
    assert(array_ptr);
    return *array_ptr;
  }

  Array &asArray()
  {
    assert(type == Type::array_ptr);
    assert(array_ptr);
    return *array_ptr;
  }

  friend bool operator==(const Value &a, const Value &b)
  {
    if (a.type != b.type) {
      cerr << "a.type: " << a.type << "\n";
      cerr << "b.type: " << b.type << "\n";
      assert(false);
    }

    switch (a.type) {
      case Type::number:
        return a.number == b.number;
      case Type::character:
        return a.character == b.character;
      case Type::array_ptr:
        return *a.array_ptr == *b.array_ptr;
    }
    assert(false);
  }

  void destroy()
  {
    switch (type) {
      case Type::number:
        destroyObject(number);
        break;
      case Type::character:
        destroyObject(character);
        break;
      case Type::array_ptr:
        destroyObject(array_ptr);
        break;
    }
  }

  void createFrom(Value &&arg)
  {
    assert(type == arg.type);

    switch (type) {
      case Type::number:
        new (&number) auto(arg.number);
        break;
      case Type::character:
        new (&character) auto(arg.character);
        break;
      case Type::array_ptr:
        assert(arg.array_ptr);
        new (&array_ptr) auto(std::move(arg.array_ptr));
        break;
    }
  }

  void createFrom(const Value &arg)
  {
    assert(type == arg.type);

    switch (type) {
      case Type::number:
        new (&number) auto(arg.number);
        break;
      case Type::character:
        new (&character) auto(arg.character);
        break;
      case Type::array_ptr:
        assert(arg.array_ptr);
        new (&array_ptr) auto(std::make_unique<Array>(*arg.array_ptr));
        break;
    }
  }

  Value& operator=(Value arg)
  {
    destroy();
    type = arg.type;
    createFrom(std::move(arg));
    return *this;
  }

  ~Value() { destroy(); }
};
}


using Values = vector<Value>;


namespace {
struct Array {
  vector<int> shape;
  Values values;

  Array() = default;
  explicit Array(const Array &) = default;
  Array(Array &&) = default;

};
}


Value::Value(Array arg)
: type(Type::array_ptr), array_ptr(std::make_unique<Array>(std::move(arg)))
{
}


static ostream& operator<<(ostream& stream, const Value &v);
static ostream& operator<<(ostream& stream, const Array &a);


static ostream& operator<<(ostream& stream, const Value &v)
{
  if (v.isNumber()) {
    stream << v.asNumber();
  }
  else if (v.isCharacter()) {
    stream << "'" << v.asCharacter() << "'";
  }
  else if (v.isArray()) {
    stream << v.asArray();
  }
  else {
    assert(false);
  }

  return stream;
}


#if 1
static ostream& operator<<(ostream& stream, const Array &a)
{
  if (a.shape.empty()) {
    stream << a.values[0];
  }
  else if (a.shape.size() == 1) {
    stream << a.values;
  }
  else if (a.shape.size() == 2) {
    stream << "[";

    for (int i=0; i!=a.shape[0]; ++i) {
      stream << " [";

      for (int j=0; j!=a.shape[1]; ++j) {
        stream << " " << a.values[i*a.shape[1] + j];
      }

      stream << " ]";
    }

    stream << " ]";
  }
  else {
    stream << "shape: " << a.shape << "\n";
    assert(false);
  }

  return stream;
}
#endif


static bool operator==(const Array &a, const Array &b)
{
  if (a.shape != b.shape) {
    if (a.values.size() == 1 && b.values.size() == 1) {
      return a.values[0] == b.values[0];
    }

    cerr << "a: " << a << "\n";
    cerr << "b: " << b << "\n";
    assert(false);
  }

  if (a.values != b.values) {
    return false;
  }

  return true;
}


namespace {
template <typename Function, typename Operator>
struct BoundOperator {
};
}


static Array makeScalarArray(Value value)
{
  if (value.isArray()) {
    return std::move(value.asArray());
  }

  Array result;
  result.values.push_back(std::move(value));
  return result;
}


static Array makeCharArray(const char *arg)
{
  Array result;
  int n = strlen(arg);
  result.shape = vector<int>{ n };
  result.values.resize(n);

  for (int i=0; i!=n; ++i) {
    result.values[i] = arg[i];
  }

  assert(result.values[0].isCharacter());
  return result;
}


namespace {
template <typename T>
struct Keyword {
};
}


namespace {
template <typename T>
struct Function {
};
}


namespace {
struct Var {
  Value *const ptr;
};
}


namespace {
template <typename T>
struct Operator {
};
}


template <typename Left, typename Right>
struct Atop {
  Left left;
  Right right;
};


template <typename Left, typename Mid, typename Right>
struct Fork {
  Left left;
  Mid mid;
  Right right;
};


namespace {
struct Shape;
struct Reshape;
struct First;
struct Each;
struct Equal;
struct Plus;
struct Times;
struct Iota;
struct Reduce;
struct Product;
struct Outer;
struct Roll;
struct Replicate;
struct MemberOf;
struct Not;
struct Empty;
struct Assign;
struct Drop;
}


template <typename Function>
static Array evaluateBinary(Array left, Array right, Function f)
{
  if (left.values.size() == 1) {
    if (right.values.size() == 1) {
      Array result;
      result.shape = vector<int>{3};

      result.values.push_back(
        f(std::move(left.values[0]), std::move(right.values[0]))
      );

      return result;
    }
    assert(false);
  }

  if (left.shape == right.shape) {
    Array result;
    result.shape = left.shape;
    result.values.resize(left.values.size());
    auto in1 = left.values.begin();
    auto in2 = right.values.begin();
    auto out = result.values.begin();

    while (out != result.values.end()) {
      *out++ = f(std::move(*in1++), std::move(*in2++));
    }

    return result;
  }

  cerr << "left.values: " << left.values << "\n";
  assert(false);
}


static Array makeArrayFromValues(Values v)
{
  if (v.empty()) {
    Array result;
    result.shape = {0};
    return result;
  }

  if (v.size() == 1) {
    return makeScalarArray(std::move(v[0]));
  }

  Array left;
  left.shape = { int(v.size()) };
  left.values = std::move(v);
  return left;
}


namespace {
struct Context {
  std::mt19937 random_engine;
};
}


static int roll(int value, Context &context)
{
  return std::uniform_int_distribution<int>(1,value)(context.random_engine);
}


namespace{
template <typename T>
struct Optional {
  bool has_value;

  Optional(T arg)
  : has_value(true),
    value(std::move(arg))
  {
  }

  Optional()
  : has_value(false)
  {
  }

  ~Optional()
  {
    if (has_value) {
      value.~T();
    }
  }

  bool operator!() const
  {
    return !has_value;
  }

  T &operator*()
  {
    assert(has_value);
    return value;
  }

  union {
    T value;
  };

  friend bool operator!=(const Optional &a, const Optional &b)
  {
    if (a.has_value != b.has_value) {
      assert(false);
    }

    if (a.has_value) {
      assert(b.has_value);
      return a.value != b.value;
    }

    assert(false);
  }
};
}


static Optional<int> maybeInteger(const Value &v)
{
  if (!v.isNumber()) {
    assert(false);
  }

  Number n = v.asNumber();
  int int_n = int(n);

  if (int_n != n) {
    assert(false);
  }

  return int_n;
}


static int product(const vector<int> &arg)
{
  return std::accumulate(arg.begin(), arg.end(), 1, std::multiplies<>());
}


static void addNCopiesTo(vector<Value> &result_values, int n, const Value &v)
{
  assert(n >= 0);

  for (int i=0; i!=n; ++i) {
    result_values.push_back(Value(v));
  }
}


static void
replicateInto(
  Values &result_values,
  const Values &left_values,
  const Values &right_values,
  int n
)
{
  for (int i=0; i!=n; ++i) {
    Optional<int> maybe_count = maybeInteger(left_values[i]);

    if (!maybe_count) {
      assert(false);
    }

    addNCopiesTo(result_values, *maybe_count, right_values[i]);
  }
}


template <typename T, typename G>
static Array evaluateNumberBinary(Fork<Array,T,Array> arg, const G &g)
{
  auto f = [&](const Value& a, const Value& b)
  {
    if (a.isNumber() && b.isNumber()) {
      return Value(g(a.asNumber(), b.asNumber()));
    }
    else {
      assert(false);
    }

    return Value();
  };

  return evaluateBinary(std::move(arg.left), std::move(arg.right), f);
}


static Value evaluate(int arg, Context &)
{
  return Number(arg);
}


static Var evaluate(Var arg, Context &)
{
  assert(arg.ptr);
  return arg;
}


static inline Value evaluate(Number arg, Context &)
{
  return arg;
}


static inline Value evaluate(char arg, Context &)
{
  return arg;
}


static Array evaluate(Value value, Context &)
{
  return makeScalarArray(std::move(value));
}


static Array evaluate(Values arg, Context &)
{
  return makeArrayFromValues(std::move(arg));
}


static Array evaluate(Atop<Function<Shape>,Array> arg, Context &)
{
  Array &array = arg.right;
  Array result;
  result.shape = { int(array.shape.size()) };

  for (int x : array.shape) {
    result.values.push_back(Number(x));
  }

  return result;
}


static Array evaluate(Atop<Function<Iota>, Array> arg, Context &)
{
  if (arg.right.shape.size() == 0) {
    int n = arg.right.values[0].asNumber();
    assert(n >= 0);
    Array result;
    result.shape = {n};

    for (int i=1; i<=n; ++i) {
      result.values.push_back(Number(i));
    }

    return result;
  }

  cerr << "arg.right: " << arg.right << "\n";
  assert(false);
}


static Array evaluate(Atop<Function<First>,Array> arg, Context &)
{
  if (arg.right.shape.size() == 0) {
    assert(false);
  }

  if (arg.right.shape.size() == 1) {
    return makeScalarArray(std::move(arg.right.values[0]));
  }

  assert(false);
}


static Array evaluate(Atop<Function<Roll>,Array> arg, Context &context)
{
  if (arg.right.shape.empty()) {
    if (arg.right.values[0].isNumber()) {
      int value = int(arg.right.values[0].asNumber());

      if (value != arg.right.values[0].asNumber()) {
        assert(false);
      }

      if (value > 0) {
        return makeScalarArray(Number(roll(value, context)));
      }
      assert(false);
    }
    assert(false);
  }

  if (arg.right.shape.size() == 1) {
    vector<Value> result;

    for (auto &x : arg.right.values) {
      if (!x.isNumber()) {
        assert(false);
      }

      int value = x.asNumber();

      if (value != arg.right.values[0].asNumber()) {
        assert(false);
      }

      result.push_back(Number(roll(value, context)));
    }

    return makeArrayFromValues(std::move(result));
  }

  cerr << "arg.right: " << arg.right << "\n";
  assert(false);
}


template <typename A, typename B, typename C>
static Atop<BoundOperator<A,B>,Function<C>>
evaluate(Atop<BoundOperator<A,B>,Function<C>> arg, Context &)
{
  return arg;
}


template <typename T>
static Function<T> evaluate(Function<T> arg, Context &)
{
  return std::move(arg);
}


template <typename T>
static Operator<T> evaluate(Operator<T> arg, Context &)
{
  return std::move(arg);
}


template <typename T>
static Keyword<T> evaluate(Keyword<T> arg, Context &)
{
  return std::move(arg);
}


static Array evaluate(Keyword<Empty>, Context &)
{
  return makeArrayFromValues({});
}


static Array evaluate(Fork<Array,Function<Equal>,Array> arg, Context &)
{
  auto f = [](auto a, auto b){ return Number(a == b); };
  return evaluateBinary(std::move(arg.left), std::move(arg.right), f);
}


static Array
evaluate(Fork<Array,Function<Drop>,Array> arg, Context &/*context*/)
{
  if (arg.left.shape.size() == 0) {
    if (arg.right.shape.size() == 1) {
      Optional<int> maybe_n_to_drop = maybeInteger(arg.left.values[0]);

      if (!maybe_n_to_drop) {
        assert(false);
      }

      int n_to_drop = *maybe_n_to_drop;

      Array result;

      if (arg.right.shape[0] <= n_to_drop) {
        result.shape = {0};
        return result;
      }

      result.shape = {arg.right.shape[0] - n_to_drop};
      auto iter = arg.right.values.begin();
      iter += n_to_drop;

      for (; iter != arg.right.values.end(); ++iter) {
        result.values.push_back(std::move(*iter));
      }

      return result;
    }
  }

  assert(false);
}


static Array evaluate(Fork<Array,Function<Reshape>,Array> arg, Context &)
{
  if (arg.left.shape.size() != 1) {
    assert(false);
  }

  Array result;

  for (Value &x : arg.left.values) {
    Optional<int> maybe_int_x = maybeInteger(x);

    if (!maybe_int_x) {
      assert(false);
    }

    result.shape.push_back(*maybe_int_x);
  }

  int n = product(result.shape);
  result.values.resize(n);

  for (int i=0; i!=n; ++i) {
    result.values[i] = Value(arg.right.values[i % arg.right.values.size()]);
  }

  return result;
}


static Array
evaluate(
  Fork<
    Array,
    Function<Fork<Operator<Outer>, Operator<Product>, Function<Times> > >,
    Array
  > arg,
  Context&
)
{
  if (arg.left.shape.size() != 1) {
    cerr << "arg.left: " << arg.left << "\n";
    cerr << "arg.left.shape.size(): " << arg.left.shape.size() << "\n";
    assert(false);
  }

  int n_rows = arg.left.shape[0];

  if (arg.right.shape.size() != 1) {
    assert(false);
  }

  int n_cols = arg.right.shape[0];

  Array result;
  result.shape = { n_rows, n_cols };
  result.values.resize(n_rows * n_cols);

  for (auto &x : arg.left.values) {
    if (!x.isNumber()) {
      assert(false);
    }
  }

  for (auto &x : arg.right.values) {
    if (!x.isNumber()) {
      assert(false);
    }
  }

  for (int i=0; i!=n_rows; ++i) {
    for (int j=0; j!=n_cols; ++j) {
      Number a = arg.left.values[i].asNumber();
      Number b = arg.right.values[j].asNumber();
      result.values[i*n_cols + j] = a * b;
    }
  }

  return result;
}


static bool elementOf(const Value &a, const Values &b)
{
  for (const Value &x : b) {
    if (a == x) {
      return true;
    }
  }

  return false;
}


static Array evaluate(Fork<Array,Function<MemberOf>,Array> arg, Context &)
{
  Array result;
  result.shape = arg.left.shape;
  size_t n = arg.left.values.size();
  result.values.resize(n);

  for (size_t i=0; i!=n; ++i) {
    result.values[i] = elementOf(arg.left.values[i], arg.right.values);
  }

  return result;
}


template <typename T>
static Array evaluate(Fork<Values,T,Array> arg, Context &context)
{
  Array left = makeArrayFromValues(std::move(arg.left));
  using Return = Fork<Array, T, Array>;
  auto x = Return{std::move(left), arg.mid, std::move(arg.right)};
  return evaluate(std::move(x), context);
}


static Array evaluate(Fork<Array,Function<Plus>,Array> arg, Context &)
{
  return
    evaluateNumberBinary(
      std::move(arg),
      [](Number a, Number b) { return a + b; }
    );
}


static Array evaluate(Fork<Array,Function<Times>,Array> arg, Context &)
{
  return
    evaluateNumberBinary(
      std::move(arg),
      [](Number a, Number b) { return a * b; }
    );
}


static Array evaluate(Array array, Context &)
{
  return array;
}


static Array evaluate(Atop<BoundOperator<Plus,Reduce>,Array> arg, Context &)
{
  const Array &right = arg.right;

  if (right.shape.size() == 1) {
    Number total = 0;

    for (auto &x : right.values) {
      if (!x.isNumber()) {
        assert(false);
      }

      total += x.asNumber();
    }

    return makeScalarArray(total);
  }

  assert(false);
}


static Array evaluate(Atop<BoundOperator<First,Each>,Array> arg, Context &)
{
  if (arg.right.shape.empty()) {
    assert(false);
  }
  else if (arg.right.shape.size() == 1) {
    Values result;

    for (auto &x : arg.right.values) {
      if (x.isNumber()) {
        result.push_back(std::move(x));
      }
      else if (x.isArray()) {
        if (x.asArray().shape.empty()) {
          assert(false);
        }
        else if (x.asArray().shape.size() == 1) {
          if (x.asArray().values.empty()) {
            assert(false);
          }

          result.push_back(std::move(x.asArray().values[0]));
        }
        else {
          assert(false);
        }
      }
      else {
        assert(false);
      }
    }

    return makeArrayFromValues(std::move(result));
  }
  else {
    assert(false);
  }

  assert(false);
}


static Array evaluate(Fork<Array, Function<Replicate>, Array> arg, Context &)
{
  Array left = std::move(arg.left);
  Array right = std::move(arg.right);
  Array result;

  if (left.shape.empty() && right.shape.empty()) {

    if (!left.values[0].isNumber()) {
      assert(false);
    }

    int n = 1;
    replicateInto(result.values, left.values, right.values, n);
  }
  else if (left.shape.size() == 1 && right.shape.size() == 1) {
    if (left.shape[0] != right.shape[0]) {
      assert(false);
    }

    int n = left.shape[0];
    replicateInto(result.values, left.values, right.values, n);
  }
  else {
    assert(false);
  }

  result.shape = { int(result.values.size()) };
  return result;
}


static Array
evaluate(
  Fork<
    Values,
    Fork<Function<Plus>,Operator<Product>,Function<Times>>,
    Array
  > arg
  , Context &
)
{
  Array left = makeArrayFromValues(std::move(arg.left));
  Array right = std::move(arg.right);

  if (left.shape.size() == 1 && right.shape.size() == 1) {
    if (left.values.size() == right.values.size()) {
      Number result = 0;
      size_t n = left.values.size();

      for (size_t i=0; i!=n; ++i) {
        assert(left.values[i].isNumber());
        assert(right.values[i].isNumber());
        result += left.values[i].asNumber()*right.values[i].asNumber();
      }

      return makeScalarArray(result);
    }
    assert(false);
  }
  assert(false);
}


static Array evaluate(Atop<Function<Not>, Array> arg, Context &)
{
  Array result;
  result.shape = arg.right.shape;
  size_t n = arg.right.values.size();
  result.values.resize(n);

  for (size_t i=0; i!=n; ++i) {
    Optional<int> maybe_x = maybeInteger(arg.right.values[0]);

    if (maybe_x!=0 && maybe_x!=1) {
      assert(false);
    }

    result.values[i] = !*maybe_x;
  }

  return result;
}


template <typename T>
static Atop<Function<T>,Array> join(Function<T> left, Array right, Context &)
{
  return { left, std::move(right) };
}


template <typename T>
static Atop<Function<T>,Array> join(Function<T> left, Value right, Context &)
{
  return { left, makeScalarArray(std::move(right)) };
}


template <typename T>
static auto join(T left, Var right, Context &context)
{
  assert(right.ptr);
  return join(left, Value(*right.ptr), context);
}


template <typename T>
static Atop<Function<T>,Array>
join(Function<T> left, Values right, Context &context)
{
  return { left, evaluate(std::move(right), context) };
}


static Value join(Var left, Atop<Keyword<Assign>, Array> right, Context&)
{
  assert(left.ptr);
  *left.ptr = std::move(right.right);
  return Value(*left.ptr);
}


static Atop<Keyword<Assign>,Array>
join(
  Keyword<Assign> left,
  Fork<Values, Function<Drop>, Array> right,
  Context& context
)
{
  return { std::move(left), evaluate(std::move(right), context) };
}


template <typename T>
static Fork<Values,T,Array> join(Value left, Atop<T,Array> right, Context &)
{
  Values new_left;
  new_left.push_back(std::move(left));
  return { std::move(new_left), std::move(right.left), std::move(right.right) };
}


template <typename T>
static Fork<Values,T,Array>
join(Value left, Fork<Values,T,Array> right, Context &)
{
  Fork<Values,T,Array> result = std::move(right);
  result.left.insert(result.left.begin(), std::move(left));
  return result;
}


template <typename T1, typename T2>
static Atop<Function<T1>,Array>
join(Function<T1> left, Atop<Function<T2>,Array> right, Context &context)
{
  return join(left, evaluate(std::move(right), context), context);
}


template <typename T>
static Atop<Operator<T>,Array>
join(Operator<T> left, Atop<Function<Iota>,Array> right, Context &context)
{
  return { left, evaluate(std::move(right), context) };
}


static Atop<Atop<Operator<Product>,Function<Times>>, Array>
join(Operator<Product>, Atop<Function<Times>,Array> right, Context &)
{
  return {{},std::move(right.right)};
}


static Atop<Fork<Function<Plus>, Operator<Product>, Function<Times>>,Array>
join(
  Function<Plus>,
  Atop<Atop<Operator<Product>, Function<Times> >, Array> right,
  Context &
)
{
  return {{},std::move(right.right)};
}


static
Atop<
  Function<Fork<Operator<Outer>, Operator<Product>, Function<Times>>>,
  Array
>
join(
  Operator<Outer> /*left*/,
  Atop<Atop<Operator<Product>, Function<Times> >, Array> right,
  Context&
)
{
  return {
    {},
    std::move(right.right)
  };
}


template <typename T>
static Atop<Operator<T>,Array>
join(Operator<T> left, Values right, Context &)
{
  return { left, makeArrayFromValues(std::move(right)) };
}


template <typename T>
static Atop<Operator<T>,Function<Iota>>
join(Operator<T>, Function<Iota>, Context &)
{
  return {};
}


static Atop<BoundOperator<Plus,Reduce>,Function<Iota>>
join(Function<Plus>, Atop<Operator<Reduce>, Function<Iota>>, Context &)
{
  return {};
}


template <typename T>
static Atop<BoundOperator<Plus,T>,Array>
join(Function<Plus> /*left*/, Atop<Operator<T>,Array> right, Context &)
{
  return { {}, std::move(right.right) };
}


template <typename T>
static Atop<BoundOperator<First,T>,Array>
join(Function<First> /*left*/, Atop<Operator<T>,Array> right, Context &)
{
  return {{},std::move(right.right)};
}


template<typename T, typename U>
static Atop<Function<T>,Array>
join(
  Function<T> left,
  Fork<Values,Function<U>,Array> right,
  Context &context
)
{
  return { std::move(left), evaluate(std::move(right), context) };
}


static Values join(Value left, Value right, Context &)
{
  Values result;
  result.push_back(std::move(left));
  result.push_back(std::move(right));
  return result;
}


static Values join(Value left, Values right, Context &)
{
  Values result;
  result.push_back(std::move(left));

  for (auto &x : right) {
    result.push_back(std::move(x));
  }

  return result;
}


static Value combine(Context &, int arg)
{
  return Value(arg);
}


template <typename T>
static Function<T> combine(Context &, Function<T> arg)
{
  return arg;
}


template <typename T>
static Keyword<T> combine(Context &, Keyword<T> arg)
{
  return arg;
}


static Array combine(Context &, Array arg)
{
  return arg;
}


static Var combine(Context &, Var arg)
{
  return arg;
}


template <typename Arg1, typename Arg2, typename ...Args>
static auto combine(Context &context, Arg1 arg1, Arg2 arg2, Args ...args)
{
  return
    join(
      evaluate(std::move(arg1), context),
      combine(context, std::move(arg2), std::move(args)...),
      context
    );
}


static Atop<BoundOperator<Plus,Reduce>,Array>
join(
  Atop<BoundOperator<Plus, Reduce>, Function<Iota> > left,
  Value right,
  Context &context
)
{
  return {
    left.left,
    evaluate(
      Atop<Function<Iota>,Array>
      {
        left.right, makeScalarArray(std::move(right))
      }
      , context
    )
  };
}


template <typename T>
struct MakeRValue {
  T operator()(T &arg) const
  {
    return T(arg);
  }
};


template <typename T>
struct MakeRValue<const T&> {
  T operator()(const T &arg) const
  {
    return T(arg);
  }
};


template <typename T>
struct MakeRValue<T&> {
  Var operator()(T &arg) const
  {
    return Var{&arg};
  }
};


template <size_t n>
struct MakeRValue<const char (&)[n]> {
  Array operator()(const char *arg) const
  {
    return makeCharArray(arg);
  }
};


namespace {
struct Placeholder {
  Context context;

  template <typename ...Args>
  auto operator()(Args &&...args)
  {
    return
      evaluate(
        combine(context, MakeRValue<Args>()(args)...),
        context
      );
  }

  static constexpr Function<Shape>     shape = {};
  static constexpr Function<Reshape>   reshape = {};
  static constexpr Function<First>     first = {};
  static constexpr Function<Equal>     equal = {};
  static constexpr Function<Plus>      plus = {};
  static constexpr Function<Times>     times = {};
  static constexpr Function<Iota>      iota = {};
  static constexpr Function<Roll>      roll = {};
  static constexpr Function<Replicate> replicate = {};
  static constexpr Function<MemberOf>  member_of = {};
  static constexpr Function<Not>       isnot = {};
  static constexpr Keyword<Empty>      empty = {};
  static constexpr Keyword<Assign>     assign = {};
  static constexpr Function<Drop>      drop = {};
  static constexpr Operator<Each>      each = {};
  static constexpr Operator<Reduce>    reduce = {};
  static constexpr Operator<Product>   product = {};
  static constexpr Operator<Outer>     outer = {};
};
}


constexpr Function<Shape>     Placeholder::shape;
constexpr Function<Reshape>   Placeholder::reshape;
constexpr Function<First>     Placeholder::first;
constexpr Function<Equal>     Placeholder::equal;
constexpr Function<Plus>      Placeholder::plus;
constexpr Function<Times>     Placeholder::times;
constexpr Function<Iota>      Placeholder::iota;
constexpr Function<Roll>      Placeholder::roll;
constexpr Function<Replicate> Placeholder::replicate;
constexpr Keyword<Empty>      Placeholder::empty;
constexpr Keyword<Assign>     Placeholder::assign;
constexpr Function<Drop>      Placeholder::drop;
constexpr Function<MemberOf>  Placeholder::member_of;
constexpr Function<Not>       Placeholder::isnot;
constexpr Operator<Each>      Placeholder::each;
constexpr Operator<Reduce>    Placeholder::reduce;
constexpr Operator<Product>   Placeholder::product;
constexpr Operator<Outer>     Placeholder::outer;


template <typename T>
static std::string str(const T &x)
{
  std::ostringstream stream;
  stream << x;
  return stream.str();
}


int main()
{
  Placeholder _;
  assert(_(1) == _(1));
  assert(!(_(2) == _(3)));
  assert(_(1,2) == _(1,2));
  assert(_(1,2,3) == _(1,2,3));
  assert(_(1,2,3).shape == vector<int>{3});
  assert(_(_.shape, 1,2,3) == _(3));
  assert(_(3, _.equal, 3) == _(1));
  assert(_(3, _.equal, _.shape, 1,2,3) == _(1));
  assert(_(_.shape,1,2,3).shape == vector<int>{1});
  assert(_(1,2,3, _.equal, 1,2,3) == _(1,1,1));
  assert(_(2, _.plus, 2) == _(4));
  assert(_(_.shape, "hello") == _(5));
  assert(_(_.iota, 3) == _(1,2,3));
  assert(_(_.plus, _.reduce, _.iota, 3) == _(6));
  assert(_(_(_.plus, _.reduce, _.iota), 3) == _(6));
  assert(_(2, _.times, 3, _.plus, 1) == _(2*(3+1)));
  assert(_(_.first, 1,2,3) == _(1));
  assert(_(_.shape, _(1,2,3), _(4,5,6)) == _(2));
  assert(_(_.shape, 1, _(2,3)) == _(2));
  assert(_(_.first, _.each, 1, 2, 3, "ABC", _(9, 8, 7)) == _(1,2,3,'A',9));
  assert(_(1,2,3,_.plus,_.product,_.times,4,5,6) == _(32));
  assert(_(_.shape, _.shape, _.roll, 6) == _(0));
  assert(_(_.shape, _.roll, 6, 6) == _(2));
  assert(_(3, _.replicate, 4) == (_(4,4,4)));
  assert(_(1, _.drop, _.iota, 5) == (_(2,3,4,5)));
  assert(_(2, _.drop, _.iota, 5) == (_(3,4,5)));
  assert(_(2, _.drop, 4, 5) == _(_.empty));
  assert(_(3, _.drop, 4, 5) == _(_.empty));
  assert(str(_(2,2, _.reshape, 1,2,3,4)) == "[ [ 1 2 ] [ 3 4 ] ]");

  assert(
    _(_(_.iota, 2), _.outer, _.product, _.times, _.iota, 2) ==
    _(2,2, _.reshape, 1, 2, 2, 4)
  );

  assert(_(1,0,2, _.replicate, 1,2,3) == _(1,3,3));

  {
    Value R = 5;
    _(R, _.assign, 1, _.drop, _.iota, 5);
    assert(*_(R).ptr == _(2,3,4,5));
  }

  assert(_(1, _.member_of, 1) == _(1));
  assert(_(1, _.member_of, 1,2,3) == _(1));
  assert(_(1,2,3, _.member_of, 2) == _(0,1,0));

  assert(_(_.isnot, 1) == _(0));

  {
    Value R = 3;
    assert(_(_.iota, R) == _(1,2,3));
  }

#if 0
  {
    Value R = 5;

    assert(
      _(
        _(_.isnot, R, _.member_of, R, _.outer, _.product, _.times, R),
        _.replicate, &R, _.assign, 1, _.drop, _.iota, R
      ) == _(2,3,5)
    );
  }
#endif
}
