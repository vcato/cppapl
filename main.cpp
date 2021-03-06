#include <cstring>
#include <random>
#include <sstream>
#include <algorithm>
#include "optional.hpp"
#include "vectorio.hpp"

using std::cerr;
using std::ostream;
using Number = double;


static int product(const vector<int> &arg)
{
  return std::accumulate(arg.begin(), arg.end(), 1, std::multiplies<>());
}


namespace {
class Array {
public:
  enum class Type { number, character, values , box };
  struct Boxed { Array &&array; };
  using Shape = vector<int>;
  using Values = vector<Array>;

public:
  Array(Number arg)
  : _type(Type::number), number(arg)
  {
  }

  Array(int arg)
  : _type(Type::number), number(arg)
  {
  }

  Array(char arg)
  : _type(Type::character), character(arg)
  {
  }

  Array(Shape shape, Values values)
  : _type(Type::values),
    box{std::move(shape), std::move(values)}
  {
    assert(!this->shape().empty());
  }

  Array(Boxed arg)
  : _type(Type::box),
    box{{}, {std::move(arg.array)}}
  {
  }

  explicit Array(const Array &arg)
  : _type(arg._type)
  {
    createFrom(arg);
  }

  Array(Array &&arg)
  : _type(arg._type)
  {
    createFrom(std::move(arg));
  }

  ~Array() { destroy(); }

  Type type() const { return _type; }

  const Shape &shape() const
  {
    assert(isNonScalar());
    return box.shape;
  }

  bool isNumber() const { return _type == Type::number; }
  bool isCharacter() const { return _type == Type::character; }
  bool isBox() const { return _type == Type::box; }
  bool isNonScalar() const { return _type == Type::values; }
  bool isSimple() const { return isNumber() || isCharacter(); }

  Number asNumber() const
  {
    assert(_type == Type::number);
    return number;
  }

  char asCharacter() const
  {
    assert(_type == Type::character);
    return character;
  }

  Values &values() &
  {
    assert(isNonScalar());
    return box.values;
  }

  const Values &values() const &
  {
    assert(isNonScalar());
    return box.values;
  }

  Array boxedValue() &&
  {
    assert(_type == Type::box);
    assert(box.shape.empty());
    return std::move(box.values[0]);
  }

  const Array& boxedValue() const &
  {
    assert(_type == Type::box);
    assert(box.shape.empty());
    return box.values[0];
  }

  Array& operator=(Array &&arg)
  {
    destroy();
    _type = arg._type;
    createFrom(std::move(arg));
    return *this;
  }

  static bool areEqual(const Array &a, const Array &b)
  {
    if (a.isNonScalar() && b.isNumber()) {
      if (a.shape().size() == 1 && a.shape()[0] == 1) {
        return areEqual(a.values()[0], b);
      }

      assert(false);
    }

    if (a._type != b._type) {
      cerr << "a.type: " << a._type << "\n";
      cerr << "b.type: " << b._type << "\n";
      assert(false);
    }

    switch (a._type) {
      case Type::number:
        return a.number == b.number;
      case Type::character:
        return a.character == b.character;
      case Type::values:
      case Type::box:
        return a.box == b.box;
    }
    assert(false);
  }

  friend ostream& operator<<(ostream &stream, Type type)
  {
    switch (type) {
      case Type::number:
        stream << "number";
        break;
      case Type::character:
        assert(false);
      case Type::values:
        stream << "non-scalar";
        break;
      case Type::box:
        stream << "box";
        break;
    }

    return stream;
  }

private:
  struct Box {
    Shape shape;
    Values values;

    Box(Shape shape, Values values)
    : shape(std::move(shape)),
      values(std::move(values))
    {
      assert(product(this->shape) == int(this->values.size()));
    }

    friend bool operator==(const Box &a, const Box &b)
    {
      return a.shape == b.shape && a.values == b.values;
    }
  };

  Type _type;

  union {
    Number number;
    char character;
    Box box;
  };

  template <typename T>
  void createFrom(T &&arg)
  {
    assert(_type == arg._type);

    switch (_type) {
      case Type::number:
        createObject(number, std::forward<T>(arg).number);
        break;
      case Type::character:
        createObject(character, std::forward<T>(arg).character);
        break;
      case Type::values:
      case Type::box:
        createObject(box, std::move(std::forward<T>(arg).box));
        break;
    }
  }

  void destroy()
  {
    switch (_type) {
      case Type::number:
        destroyObject(number);
        break;
      case Type::character:
        destroyObject(character);
        break;
      case Type::values:
      case Type::box:
        destroyObject(box);
        break;
    }
  }
};
}


using Values = vector<Array>;


namespace {
static bool operator==(const Array &a, const Array &b)
{
  return Array::areEqual(a,b);
}
}


static ostream& operator<<(ostream& stream, const Array &a);


static bool isAllCharacters(const Values &v, int start, int n)
{
  for (auto it = v.begin() + start; it != v.begin() + start + n; ++it) {
    const Array &a = *it;

    if (!a.isCharacter()) {
      return false;
    }
  }

  return true;
}


static bool isAllCharacters(const Values &v)
{
  return isAllCharacters(v, 0, v.size());
}


static void
printSpanOn(ostream &stream, const Values &v, int start, int n)
{
  if (isAllCharacters(v, start, n)) {
    stream << " \"";

    for (int i=start; i!=start+n; ++i) {
      stream << v[i].asCharacter();
    }

    stream << "\"";
  }
  else {
    stream << " [";

    for (int i=start; i!=start+n; ++i) {
      stream << " " << v[i];
    }

    stream << " ]";
  }
}


static void
printNonScalarArrayOn(
  ostream &stream, const Array::Shape &shape, const Array::Values &values
)
{
  assert(!shape.empty());

  if (shape.size() == 1 && isAllCharacters(values)) {
    printSpanOn(stream, values, 0, values.size());
    return;
  }

  stream << "(";

  for (int i=0; i!=shape[0]; ++i) {
    if (shape.size() == 1) {
      stream << " " << values[i];
    }
    else if (shape.size() == 2) {
      printSpanOn(stream, values, i*shape[1], shape[1]);
    }
    else {
      stream << "shape: " << shape << "\n";
      assert(false);
    }
  }

  stream << " )";
}


static ostream& operator<<(ostream& stream, const Array &v)
{
  if (v.isNumber()) {
    stream << v.asNumber();
  }
  else if (v.isCharacter()) {
    stream << "'" << v.asCharacter() << "'";
  }
  else if (v.isNonScalar()) {
    printNonScalarArrayOn(stream, v.shape(), v.values());
  }
  else if (v.isBox()) {
    stream << "box("  << v.boxedValue() << ")";
  }
  else {
    assert(false);
  }

  return stream;
}


static Array makeArrayFromString(const char *arg)
{
  int n = strlen(arg);
  Array::Shape shape = { n };
  Values values;

  for (int i=0; i!=n; ++i) {
    values.push_back(arg[i]);
  }

  assert(values[0].isCharacter());
  return Array(std::move(shape), std::move(values));
}


namespace {
template <typename T> struct Keyword { };
template <typename T> struct Function { };
template <typename T> struct Dfn { T f; };
template <typename T> struct Index { T arg; };
template <typename T> struct Operator { };

template <typename Function, typename Operator>
struct BoundOperator {
  Function left;
};

}



namespace {
struct Context {
  std::mt19937 &random_engine;
  Array *right_arg_ptr = nullptr;
};
}


namespace {
template <typename F>
struct Expr {
  Context &context;
  F f;

  Expr(Context &context, F f)
  : context(context), f(std::move(f))
  {
  }

  bool evaluated = false;

  Expr(const Expr &) = delete;
  void operator=(const Expr &) = delete;
  void operator=(Expr &&) = delete;

  Expr(Expr &&arg)
  : context(arg.context),
    f(std::move(arg.f))
  {
    arg.evaluated = true;
  }

  ~Expr();

  operator Array() &&;
};
}


namespace {
struct Var { Array *const ptr; };
}


using Vars = vector<Var>;


namespace {
template <typename Left, typename Right>
struct Atop {
  Left left;
  Right right;
};
}


namespace {
template <typename Left, typename Right>
struct Partial {
  Left left;
  Right right;
};
}


namespace {
template <typename Left, typename Mid, typename Right>
struct Fork {
  Left left;
  Mid mid;
  Right right;
};
}


namespace {
struct Shape;
struct Reshape;
struct First;
struct Each;
struct Equal;
struct Plus;
struct Minus;
struct Times;
struct Divide;
struct Power;
struct Greater;
struct Iota;
struct Reduce;
struct Product;
struct Outer;
struct Beside;
struct Roll;
struct Replicate;
struct Enclose;
struct GradeUp;
struct MemberOf;
struct Not;
struct Empty;
struct Assign;
struct Drop;
struct RightArg;
struct Where;
}


static bool isScalar(const Array &a)
{
  if (a.isNonScalar()) {
    assert(!a.shape().empty());
    return false;
  }
  else {
    return true;
  }
}


template <typename A, typename B>
static Atop<A,B> atop(A a, B b)
{
  return { std::move(a), std::move(b) };
}


template <typename A, typename B, typename C>
static Fork<A,B,C> fork(A a, B b, C c)
{
  return { std::move(a), std::move(b), std::move(c) };
}


template <typename Function>
static Array evaluateUnary(Array a, Function f)
{
  if (isScalar(a)) {
    return f(a);
  }

  size_t n = a.values().size();
  Values values;

  for (size_t i=0; i!=n; ++i) {
    values.push_back(f(a.values()[i]));
  }

  return Array(a.shape(), std::move(values));
}


template <typename Function>
static Array evaluateBinary(Array left, Array right, Function f)
{
  if (isScalar(left)) {
    if (isScalar(right)) {
      if (left.isSimple()) {
        if (right.isSimple()) {
          return f(std::move(left), std::move(right));
        }
        else if (right.isBox()) {
          Array result =
            evaluateBinary(std::move(left), std::move(right).boxedValue(), f);

          return result;
        }
        else {
          assert(false);
        }
      }
      else {
        assert(false);
      }
    }

    assert(!right.shape().empty());
    Values values;

    for (auto &x : right.values()) {
      values.push_back(evaluateBinary(Array(left), std::move(x), f));
    }

    return Array(right.shape(), std::move(values));
  }

  if (isScalar(right)) {
    Values values;

    for (auto &x : left.values()) {
      values.push_back(evaluateBinary(std::move(x), Array(right), f));
    }

    return Array(left.shape(), std::move(values));
  }

  if (left.values().size() == 1) {
    if (right.values().size() == 1) {
      Array::Shape shape = {1};
      Values values;

      values.push_back(
        f(std::move(left.values()[0]), std::move(right.values()[0]))
      );

      return Array(std::move(shape), std::move(values));
    }
    assert(false);
  }

  if (left.shape() == right.shape()) {
    Array::Shape shape = left.shape();
    Values values;
    size_t n = left.values().size();
    auto in1 = left.values().begin();
    auto in2 = right.values().begin();

    for (size_t i=0; i!=n; ++i) {
      values.push_back(f(std::move(*in1++), std::move(*in2++)));
    }

    return Array(std::move(shape), std::move(values));
  }

  cerr << "left.values: " << left.values() << "\n";
  assert(false);
}


static Array makeArrayFromValues(Values v)
{
  if (v.empty()) {
    return Array({0}, {});
  }

  size_t n = v.size();

  if (n == 1) {
    return std::move(v[0]);
  }

  return Array({ int(n) }, std::move(v));
}


static int roll(int range, Context &context)
{
  return std::uniform_int_distribution<int>(1,range)(context.random_engine);
}


static Optional<int> maybeInteger(const Array &a)
{
  if (isScalar(a)) {
    if (a.isNumber()) {
      Number n = a.asNumber();
      int int_n = int(n);

      if (int_n != n) {
        assert(false);
      }

      return int_n;
    }
    else {
      assert(false);
    }
  }
  else {
    return {};
  }
}


static void addNCopiesTo(Values &result_values, int n, const Array &v)
{
  assert(n >= 0);

  for (int i=0; i!=n; ++i) {
    result_values.push_back(Array(v));
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
  auto f = [&](const Array& a, const Array& b)
  {
    if (a.isNumber()) {
      if (b.isNumber()) {
        return g(a.asNumber(), b.asNumber());
      }
      else if (b.isBox()) {
        assert(false);
      }
      else {
        assert(false);
      }
    }
    else {
      assert(false);
    }
  };

  return evaluateBinary(std::move(arg.left), std::move(arg.right), f);
}


static bool isElementOf(const Array &a, const Array &b)
{
  assert(isScalar(a));

  if (isScalar(b)) {
    return (a == Array(b));
  }

  for (const Array &x : b.values()) {
    if (a == x) {
      return true;
    }
  }

  return false;
}


static bool isVector(const Array &arg)
{
  return arg.isNonScalar() && arg.shape().size() == 1;
}


static Optional<Array> maybeElement(const Array &a, const Array &index_array)
{
  if (Optional<int> maybe_index = maybeInteger(index_array)) {
    int index = *maybe_index;

    if (index >= 1 && index <= a.shape()[0]) {
      return Array(a.values()[index-1]);
    }
    else {
      assert(false);
    }
  }
  else if (isVector(index_array)) {
    if (a.shape().size() == 2) {
      if (index_array.shape()[0] == 2) {
        Optional<int> maybe_i1 = maybeInteger(index_array.values()[0]);
        Optional<int> maybe_i2 = maybeInteger(index_array.values()[1]);

        if (!maybe_i1 || !maybe_i2) {
          assert(false);
        }

        int i1 = *maybe_i1;
        int i2 = *maybe_i2;

        if (i1 < 1 || i1 > a.shape()[0]) {
          assert(false);
        }

        if (i2 < 1 || i2 > a.shape()[1]) {
          assert(false);
        }

        return Array(a.values()[(i1-1)*a.shape()[1] + (i2-1)]);
      }
      assert(false);
    }

    assert(false);
  }
  else if (index_array.isNonScalar()) {
    Values values;

    for (const auto &x : index_array.values()) {
      Optional<Array> y = maybeElement(a,x);

      if (!y) {
        return {};
      }

      values.push_back(std::move(*y));
    }

    return Array(index_array.shape(), std::move(values));
  }
  else {
    assert(false);
  }
}


static bool isNot(const Array &a)
{
  Optional<int> maybe_x = maybeInteger(a);

  if (maybe_x != 0 && maybe_x != 1) {
    assert(false);
  }

  return !*maybe_x;
}


namespace {
Array evaluate(int arg, Context &)
{
  return Number(arg);
}
}


namespace {
Var evaluate(Var arg, Context &)
{
  assert(arg.ptr);
  return arg;
}
}


namespace {
vector<Var> evaluate(vector<Var> arg, Context &)
{
  return arg;
}
}


namespace {
inline Array evaluate(Number arg, Context &)
{
  return arg;
}
}


namespace {
inline Array evaluate(char arg, Context &)
{
  return arg;
}
}


namespace {
Array evaluate(Array value, Context &)
{
  return value;
}
}


namespace {
Array evaluate(Values arg, Context &)
{
  return makeArrayFromValues(std::move(arg));
}
}


namespace {
Array evaluate(Atop<Function<Shape>,Array> arg, Context &)
{
  Values values;

  if (isScalar(arg.right)) {
    return Array({0},{});
  }

  const Array::Shape &shape = arg.right.shape();

  for (int x : shape) {
    values.push_back(Number(x));
  }

  return Array({ int(shape.size()) }, std::move(values));
}
}


namespace {
Array evaluate(Atop<Function<Divide>,Array> arg, Context &)
{
  auto f = [](const Array &a){
    if (!a.isNumber()) {
      assert(false);
      return Array(0);
    }

    return Array(1/a.asNumber());
  };

  return evaluateUnary(std::move(arg.right), f);
}
}


namespace {
Array evaluate(Atop<Function<Iota>, Array> arg, Context &)
{
  if (isScalar(arg.right)) {
    Optional<int> maybe_n = maybeInteger(arg.right);

    if (!maybe_n) {
      assert(false);
    }

    int n = *maybe_n;
    assert(n >= 0);
    Array::Shape shape = {n};

    Values values;

    for (int i=1; i<=n; ++i) {
      values.push_back(Number(i));
    }

    return Array(std::move(shape), std::move(values));
  }

  if (isVector(arg.right)) {
    if (arg.right.shape()[0] == 2) {
      Optional<int> maybe_n = maybeInteger(arg.right.values()[0]);
      Optional<int> maybe_m = maybeInteger(arg.right.values()[1]);

      if (!maybe_n) {
        assert(false);
      }

      int n = *maybe_n;

      if (!maybe_m) {
        assert(false);
      }

      int m = *maybe_m;
      Values values;

      for (int i=0; i!=n; ++i) {
        for (int j=0; j!=m; ++j) {
          values.push_back(makeArrayFromValues({i+1, j+1}));
        }
      }

      Array::Shape shape = {n,m};
      return Array(std::move(shape), std::move(values));
    }

    assert(false);
  }

  cerr << "arg.right: " << arg.right << "\n";
  assert(false);
}
}


namespace {
Array evaluate(Atop<Function<First>,Array> arg, Context &)
{
  if (isScalar(arg.right)) {
    if (arg.right.isBox()) {
      return std::move(arg.right).boxedValue();
    }
    assert(false);
  }

  if (arg.right.shape().size() == 1) {
    return std::move(arg.right.values()[0]);
  }

  assert(false);
}
}


namespace {
Array evaluate(Atop<Function<Roll>,Array> arg, Context &context)
{
  if (isScalar(arg.right)) {
    Optional<int> maybe_value = maybeInteger(arg.right);

    if (maybe_value) {
      int value = *maybe_value;

      if (value > 0) {
        return Number(roll(value, context));
      }

      assert(false);
    }

    assert(false);
  }

  if (arg.right.shape().size() == 1) {
    Values result;

    for (auto &x : arg.right.values()) {
      if (!x.isNumber()) {
        assert(false);
      }

      int value = x.asNumber();

      if (value != arg.right.values()[0].asNumber()) {
        assert(false);
      }

      result.push_back(Number(roll(value, context)));
    }

    return makeArrayFromValues(std::move(result));
  }

  cerr << "arg.right: " << arg.right << "\n";
  assert(false);
}
}


namespace {
Array evaluate(Fork<Array, Function<Roll>, Array> arg, Context &context)
{
  if (arg.left.isNumber() && arg.right.isNumber()) {
    Optional<int> maybe_count = maybeInteger(arg.left);
    Optional<int> maybe_range = maybeInteger(arg.right);

    if (maybe_count && maybe_range) {
      int count = *maybe_count;
      int range = *maybe_range;
      Values values;

      for (int i=0; i!=count; ++i) {
        values.push_back(roll(range, context));
      }

      return Array({count}, std::move(values));
    }

    assert(false);
  }

  cerr << "arg.left: " << arg.left << "\n";
  cerr << "arg.right: " << arg.right << "\n";
  assert(false);
}
}


namespace {
Array evaluate(Fork<Array,Function<Greater>,Array> arg, Context &)
{
  auto f = [](Number a, Number b){ return a > b; };
  return evaluateNumberBinary(arg, f);
}
}


namespace {
template <typename A, typename B, typename C>
Atop<BoundOperator<A,B>,Function<C>>
evaluate(Atop<BoundOperator<A,B>,Function<C>> arg, Context &)
{
  return arg;
}
}


namespace {
template <typename T>
Function<T> evaluate(Function<T> arg, Context &)
{
  return arg;
}
}


namespace {
template <typename T>
Operator<T> evaluate(Operator<T> arg, Context &)
{
  return arg;
}
}


namespace {
template <typename T>
Keyword<T> evaluate(Keyword<T> arg, Context &)
{
  return arg;
}
}


namespace {
Array evaluate(Keyword<RightArg> /*arg*/, Context &context)
{
  if (!context.right_arg_ptr) {
    assert(false);
  }

  return Array(*context.right_arg_ptr);
}
}


namespace {
Array evaluate(Keyword<Empty>, Context &)
{
  return makeArrayFromValues({});
}
}


namespace {
Array evaluate(Fork<Array,Function<Equal>,Array> arg, Context &)
{
  auto f = [](auto a, auto b){ return Number(a == b); };
  return evaluateBinary(std::move(arg.left), std::move(arg.right), f);
}
}


namespace {
Array evaluate(Fork<Array,Function<Drop>,Array> arg, Context &/*context*/)
{
  if (isScalar(arg.left)) {
    if (arg.right.shape().size() == 1) {
      Optional<int> maybe_n_to_drop = maybeInteger(arg.left);

      if (!maybe_n_to_drop) {
        assert(false);
      }

      int n_to_drop = *maybe_n_to_drop;

      if (arg.right.shape()[0] <= n_to_drop) {
        return Array({0}, {});
      }

      Array::Shape shape = {arg.right.shape()[0] - n_to_drop};
      auto iter = arg.right.values().begin();
      iter += n_to_drop;
      Values values;

      for (; iter != arg.right.values().end(); ++iter) {
        values.push_back(std::move(*iter));
      }

      return Array(std::move(shape), std::move(values));
    }
  }

  assert(false);
}
}


namespace {
Array evaluate(Fork<Array,Function<Reshape>,Array> arg, Context &)
{
  if (Optional<int> maybe_n = maybeInteger(arg.left)) {
    int n = *maybe_n;

    if (isScalar(arg.right)) {
      Values values;

      for (int i=0; i!=n; ++i) {
        values.push_back(arg.right);
      }

      Array::Shape shape = {n};
      return Array(std::move(shape), std::move(values));
    }
    assert(false);
  }

  Array::Shape shape;

  for (Array &x : arg.left.values()) {
    Optional<int> maybe_int_x = maybeInteger(x);

    if (!maybe_int_x) {
      assert(false);
    }

    shape.push_back(*maybe_int_x);
  }

  int n = product(shape);
  Values values;

  for (int i=0; i!=n; ++i) {
    values.push_back(arg.right.values()[i % arg.right.values().size()]);
  }

  return Array(std::move(shape), std::move(values));
}
}


namespace {
Array
evaluate(
  Fork<
    Array,
    Function<Fork<Operator<Outer>, Operator<Product>, Function<Times> > >,
    Array
  > arg,
  Context&
)
{
  if (arg.left.shape().size() != 1) {
    cerr << "arg.left: " << arg.left << "\n";
    cerr << "arg.left.shape.size(): " << arg.left.shape().size() << "\n";
    assert(false);
  }

  int n_rows = arg.left.shape()[0];

  if (arg.right.shape().size() != 1) {
    assert(false);
  }

  int n_cols = arg.right.shape()[0];

  Array::Shape shape = { n_rows, n_cols };
  Values values;

  for (auto &x : arg.left.values()) {
    if (!x.isNumber()) {
      assert(false);
    }
  }

  for (auto &x : arg.right.values()) {
    if (!x.isNumber()) {
      assert(false);
    }
  }

  for (int i=0; i!=n_rows; ++i) {
    for (int j=0; j!=n_cols; ++j) {
      Number a = arg.left.values()[i].asNumber();
      Number b = arg.right.values()[j].asNumber();
      values.push_back(a * b);
    }
  }

  return Array(std::move(shape), std::move(values));
}
}


namespace {
Array evaluate(Fork<Array,Function<MemberOf>,Array> arg, Context &)
{
  if (isScalar(arg.left)) {
    Array left_value = std::move(arg.left);
    return isElementOf(left_value, arg.right);
  }

  size_t n = arg.left.values().size();
  Values values;

  for (size_t i=0; i!=n; ++i) {
    values.push_back(isElementOf(arg.left.values()[i], arg.right));
  }

  return Array(arg.left.shape(), std::move(values));
}
}


namespace {
Array elements(const Array &a, const Array &indices)
{
  int n = indices.shape()[0];
  Values values;

  for (int i=0; i!=n; ++i) {
    const Array &index = indices.values()[i];

    Optional<Array> maybe_element = maybeElement(a, index);

    if (maybe_element) {
      values.push_back(std::move(*maybe_element));
    }
    else {
      assert(false);
    }
  }

  return Array({n}, std::move(values));
}
}


namespace {
Array evaluate(Atop<Array,Index<Array>> arg, Context &/*context*/)
{
  Array &right = arg.right.arg;

  if (isVector(arg.left)) {
    if (isVector(right)) {
      return elements(/*array*/arg.left, /*indices*/right);
    }

    if (Optional<Array> maybe_element = maybeElement(arg.left, right)) {
      return std::move(*maybe_element);
    }

    cerr << "arg.right: " << right << "\n";
    assert(false);
  }

  if (arg.left.shape().size() == 2) {
    if (isVector(right)) {
      return elements(/*array*/arg.left, /*indices*/right);
    }

    assert(false);
  }

  cerr << "arg.left: " << arg.left << "\n";
  cerr << "right: " << right << "\n";
  assert(false);
}
}


namespace {
template <typename T>
Array evaluate(Fork<Values,T,Array> arg, Context &context)
{
  Array left = makeArrayFromValues(std::move(arg.left));
  using Return = Fork<Array, T, Array>;
  auto x = Return{std::move(left), arg.mid, std::move(arg.right)};
  return evaluate(std::move(x), context);
}
}


namespace {
template <typename T>
Array evaluate(Atop<Values,T> arg, Context &context)
{
  return
    evaluate(
      atop(
        makeArrayFromValues(std::move(arg.left)),
        std::move(arg.right)
      ),
      context
    );
}
}


namespace {
Array evaluate(Fork<Array,Function<Plus>,Array> arg, Context &)
{
  auto f = [](Number a, Number b) { return a + b; };
  return evaluateNumberBinary(std::move(arg), f);
}
}


namespace {
Array evaluate(Fork<Array,Function<Minus>,Array> arg, Context &)
{
  auto f = [](Number a, Number b) { return a - b; };
  return evaluateNumberBinary(std::move(arg), f);
}
}


namespace {
Array evaluate(Fork<Array,Function<Times>,Array> arg, Context &)
{
  auto f = [](Number a, Number b) { return a * b; };
  return evaluateNumberBinary(std::move(arg), f);
}
}


namespace {
Array evaluate(Fork<Array, Function<Divide>, Array> arg, Context&)
{
  auto f = [](Number a, Number b) { return a / b; };
  return evaluateNumberBinary(arg, f);
}
}


namespace {
Array evaluate(Fork<Array, Function<Power>, Array> arg, Context&)
{
  auto f = [](Number a, Number b) { return std::pow(a,b); };
  return evaluateNumberBinary(arg, f);
}
}


namespace {
Array
evaluate(Atop<BoundOperator<Function<Plus>,Reduce>,Array> arg, Context &)
{
  const Array &right = arg.right;

  if (right.shape().size() == 1) {
    Number total = 0;

    for (auto &x : right.values()) {
      if (!x.isNumber()) {
        assert(false);
      }

      total += x.asNumber();
    }

    return total;
  }

  assert(false);
}
}


namespace {
template <typename T>
Array
evaluate(
  Atop<BoundOperator<Dfn<T> , Each>, Array> arg,
  Context& context
)
{
  if (isScalar(arg.right)) {
    assert(false);
  }
  else if (arg.right.shape().size() == 1) {
    assert(false);
  }
  else {
    Values values;

    for (auto &x : arg.right.values()) {
      values.push_back(evaluateDfn(arg.left.left, std::move(x), context));
    }

    return Array(std::move(arg).right.shape(), std::move(values));
  }

  assert(false);
}
}


namespace {
Array
evaluate(Atop<BoundOperator<Function<First>,Each>,Array> arg, Context &)
{
  if (arg.right.shape().empty()) {
    assert(false);
  }
  else if (arg.right.shape().size() == 1) {
    Values result;

    for (auto &x : arg.right.values()) {
      if (x.isNumber()) {
        result.push_back(std::move(x));
      }
      else if (x.isNonScalar()) {
        if (x.shape().empty()) {
          assert(false);
        }
        else if (x.shape().size() == 1) {
          if (x.values().empty()) {
            assert(false);
          }

          result.push_back(std::move(x.values()[0]));
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
}


namespace {
Array evaluate(Fork<Array, Function<Replicate>, Array> arg, Context &)
{
  Array left = std::move(arg.left);
  Array right = std::move(arg.right);
  Values values;

  if (isScalar(left) && isScalar(right)) {
    Optional<int> maybe_n = maybeInteger(left);

    if (!maybe_n) {
      assert(false);
    }

    int n = *maybe_n;

    if (n < 0) {
      assert(false);
    }

    Array::Shape shape = {n};
    Values values;

    for (int i=0; i!=n; ++i) {
      values.push_back(Array(right));
    }

    return Array(std::move(shape), std::move(values));
  }
  else if (left.shape().size() == 1 && right.shape().size() == 1) {
    if (left.shape()[0] != right.shape()[0]) {
      assert(false);
    }

    int n = left.shape()[0];
    replicateInto(values, left.values(), right.values(), n);
  }
  else {
    assert(false);
  }

  size_t n = values.size();
  return Array({ int(n) }, std::move(values));
}
}


namespace {
Array
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

  if (left.shape().size() == 1 && right.shape().size() == 1) {
    if (left.values().size() == right.values().size()) {
      Number result = 0;
      size_t n = left.values().size();

      for (size_t i=0; i!=n; ++i) {
        assert(left.values()[i].isNumber());
        assert(right.values()[i].isNumber());
        result += left.values()[i].asNumber()*right.values()[i].asNumber();
      }

      return result;
    }
    assert(false);
  }
  assert(false);
}
}


namespace {
Array evaluate(Atop<Function<Not>, Array> arg, Context &)
{
  auto f = [](const Array& x){ return isNot(x); };
  return evaluateUnary(std::move(arg.right), f);
}
}


namespace {
Array evaluate(Atop<Function<Enclose>,Array> arg, Context &)
{
  if (arg.right.isSimple()) {
    return std::move(arg.right);
  }

  return Array::Boxed{std::move(arg.right)};
}
}


namespace {
Array evaluate(Atop<Function<GradeUp>, Array> arg, Context &)
{
  if (isVector(arg.right)) {
    Array::Shape indices;
    int n = arg.right.shape()[0];

    for (int i=0; i!=n; ++i) {
      indices.push_back(i);
    }

    const Values &v = arg.right.values();

    auto comp = [&](int a, int b)
    {
      if (v[a].isNumber() && v[b].isNumber()) {
        return v[a].asNumber() < v[b].asNumber();
      }
      else {
        assert(false);
        return true;
      }
    };

    std::sort(indices.begin(), indices.end(), comp);
    Values values;

    for (auto &x : indices) {
      values.push_back(x+1);
    }

    return Array({n}, std::move(values));
  }

  cerr << "arg.right: " << arg.right << "\n";
  assert(false);
}
}


namespace {
Array
evaluate(
  Fork<
    Array,
    Fork<
      Function<Plus>,
      Operator<Beside>,
      Function<Divide>
    >,
    Array
  > arg,
  Context& context
)
{
  Array a1 = std::move(arg.left);
  Array a2 = std::move(arg.right);
  auto f1 = std::move(arg.mid.left);
  auto f2 = std::move(arg.mid.right);
  Array a = evaluate(atop(std::move(f2), std::move(a2)), context);
  Array b = evaluate(fork(std::move(a1),std::move(f1),std::move(a)), context);
  return b;
}
}


template <typename T>
struct MakeRValue {
  T operator()(T arg) const
  {
    return arg;
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
    return makeArrayFromString(arg);
  }
};


static Array combine(Context &, int arg)
{
  return Array(arg);
}


static Array combine(Context &, Number arg)
{
  return Array(arg);
}


static Array combine(Context &context, Keyword<RightArg>)
{
  if (!context.right_arg_ptr) {
    assert(false);
  }

  return Array(*context.right_arg_ptr);
}


namespace {
template <typename T>
Atop<Function<T>,Array> join(Function<T> left, Array right, Context &)
{
  return { left, std::move(right) };
}
}


namespace {
Partial<Keyword<Assign>,Array>
join(Keyword<Assign> left, Array right, Context &)
{
  return { std::move(left), std::move(right) };
}
}


namespace {
template <typename T>
Partial<Operator<T>,Array>
join(Operator<T> left, Array right, Context &)
{
  return { std::move(left), std::move(right) };
}
}


namespace {
template <typename T, typename U>
Partial<Operator<T>,Array>
join(Operator<T> left, Fork<Values, Function<U>, Array> right, Context &context)
{
  return { std::move(left), evaluate(std::move(right), context) };
}
}


namespace {
template <typename T>
auto join(T left, Var right, Context &context)
{
  assert(right.ptr);
  return join(std::move(left), Array(*right.ptr), context);
}
}


namespace {
Vars join(Var left, Var right, Context &)
{
  return { std::move(left), std::move(right) };
}
}


namespace {
template <typename T>
Atop<Function<T>,Array>
join(Function<T> left, Values right, Context &context)
{
  return { left, evaluate(std::move(right), context) };
}
}


namespace {
Fork<Var,Keyword<Assign>,Array>
join(Var left, Partial<Keyword<Assign>, Array> right, Context&)
{
  return { left, std::move(right.left), std::move(right.right) };
}
}


namespace {
Array
evaluate(
  Atop<
    Fork< vector<Array>, Operator<Beside>, Function<Reshape> >,
    Array
  > arg,
  Context& context
)
{
  return
    evaluate(
      fork(
        std::move(arg.left.left),
        std::move(arg.left.right),
        std::move(arg.right)
      ),
      context
    );
}
}


namespace {
Array
evaluate(
  Fork<
    Fork< vector<Array>, Operator<Beside>, Function<Reshape> >,
    Operator<Each>,
    Array
  > arg,
  Context& context
)
{
  if (isScalar(arg.right)) {
    assert(false);
  }

  if (arg.right.shape().size() == 1) {
    Values values;

    for (auto &x : arg.right.values()) {
      values.push_back(evaluate(atop(arg.left, std::move(x)), context));
    }

    return Array(std::move(arg.right.shape()), std::move(values));
  }

  assert(false);
}
}


namespace {
template <typename T>
Partial<Keyword<Assign>,Array>
join(Keyword<Assign> left, T right, Context& context)
{
  return { std::move(left), evaluate(std::move(right), context) };
}
}


namespace {
template <typename T>
Fork<Values,T,Array> join(Array left, Atop<T,Array> right, Context &)
{
  Values new_left;
  new_left.push_back(std::move(left));
  return { std::move(new_left), std::move(right.left), std::move(right.right) };
}
}


namespace {
template <typename T>
Fork<Values,T,Array>
join(Array left, Fork<Values,T,Array> right, Context &)
{
  Fork<Values,T,Array> result = std::move(right);
  result.left.insert(result.left.begin(), std::move(left));
  return result;
}
}


namespace {
template <typename T1, typename T2>
Atop<Function<T1>,Array>
join(Function<T1> left, Atop<Function<T2>,Array> right, Context &context)
{
  return join(left, evaluate(std::move(right), context), context);
}
}


namespace {
template <typename T>
Fork<Keyword<RightArg>, Function<T>, Array>
join(Keyword<RightArg> left, Atop<Function<T>, Array> right, Context &)
{
  return { std::move(left), std::move(right.left), std::move(right.right) };
}
}


namespace {
Atop< BoundOperator<Function<Plus>, Reduce>, Array >
join(
  Function<Plus> left,
  Atop<
    Partial<
      Operator<Reduce>,
      Function<Iota>
    >,
    Array
  > right,
  Context& context
)
{
  return {
    { std::move(left) },
    evaluate(
      atop(
        std::move(right.left.right),
        std::move(right.right)
      ),
      context
    )
  };
}
}


namespace {
Atop<Fork<Function<Plus>,Operator<Beside>,Function<Divide>>, Array>
join(
  Function<Plus> left,
  Atop<
    Partial<
      Operator<Beside>,
      Function<Divide>
    >,
    Array
  > right,
  Context&
)
{
  return {
    {
      std::move(left),
      std::move(right.left.left),
      std::move(right.left.right)
    },
    std::move(right.right)
  };
}
}


namespace {
template <typename T, typename U>
Atop<Partial<Operator<U>,Function<T>>, Array>
join(Operator<U>, Atop<Function<T>,Array> right, Context &)
{
  return {{},std::move(right.right)};
}
}


namespace {
Atop<Fork<Function<Plus>, Operator<Product>, Function<Times>>,Array>
join(
  Function<Plus>,
  Atop<Partial<Operator<Product>, Function<Times> >, Array> right,
  Context &
)
{
  return {{},std::move(right.right)};
}
}


namespace {
Atop<
  Function<Fork<Operator<Outer>, Operator<Product>, Function<Times>>>,
  Array
>
join(
  Operator<Outer> /*left*/,
  Atop<Partial<Operator<Product>, Function<Times> >, Array> right,
  Context&
)
{
  return {
    {},
    std::move(right.right)
  };
}
}


namespace {
template <typename T>
Partial<Operator<T>,Array>
join(Operator<T> left, Values right, Context &)
{
  return { left, makeArrayFromValues(std::move(right)) };
}
}


namespace {
template <typename T, typename U>
Partial<Operator<T>,Function<U>> join(Operator<T>, Function<U>, Context &)
{
  return {};
}
}


namespace {
template <typename T, typename U, typename V>
Atop<BoundOperator<Function<T>,V>,Function<U>>
join(Function<T>, Partial<Operator<V>, Function<U>>, Context &)
{
  return {};
}
}


namespace {
template<typename T, typename U>
Atop<Function<T>,Array>
join(
  Function<T> left,
  Fork<Values,Function<U>,Array> right,
  Context &context
)
{
  return { std::move(left), evaluate(std::move(right), context) };
}
}


namespace {
Values join(Array left, Array right, Context &)
{
  Values result;
  result.push_back(std::move(left));
  result.push_back(std::move(right));
  return result;
}
}


namespace {
template <typename T>
Fork<Values,T,Array> join(Array left, Partial<T,Array> right, Context &)
{
  Values new_left;
  new_left.push_back(std::move(left));
  return { std::move(new_left), std::move(right.left), std::move(right.right) };
}
}


namespace {
Values join(Array left, Values right, Context &)
{
  Values result;
  result.push_back(std::move(left));

  for (auto &x : right) {
    result.push_back(std::move(x));
  }

  return result;
}
}


namespace {
template <typename T>
auto join(Var left, T right, Context &context)
{
  return join(Array(*left.ptr), std::move(right), context);
}
}


namespace {
template <typename T>
Atop<Dfn<T>,Array>
join(Dfn<T> left, Array right, Context&)
{
  return {std::move(left), std::move(right)};
}
}


namespace {
Atop<Values, Index<Array>>
join(Array left, Atop<Array, Index<Array>> right, Context&)
{
  return {
    { std::move(left), std::move(right.left) },
    std::move(right.right)
  };
}
}


namespace {
Atop<Values, Index<Array>>
join(Array left, Atop<Values, Index<Array>> right, Context&)
{
  Values new_left;
  new_left.push_back(std::move(left));

  for (auto &x : right.left) {
    new_left.push_back(std::move(x));
  }

  return { std::move(new_left), std::move(right.right) };
}
}


static Array makeArrayFromVars(Vars vars)
{
  Values values;

  for (auto &var : vars) {
    values.push_back(Array(*var.ptr));
  }

  return makeArrayFromValues(std::move(values));
}


namespace {
template <typename T>
auto join(Function<T> left, Vars right, Context &context)
{
  return join(std::move(left), makeArrayFromVars(std::move(right)), context);
}
}


namespace {
template <typename T, typename U, typename V>
auto
join(
  Function<T> left,
  Atop<BoundOperator<U, V>, Array> right,
  Context &context
)
{
  return join(left, evaluate(std::move(right), context), context);
}
}


namespace {
template <typename T>
Atop< BoundOperator<Dfn<T>,Each> , Array >
join(
  Dfn<T> left,
  Partial<Operator<Each>, Array> right,
  Context&
)
{
  return { { std::move(left) }, std::move(right.right) };
}
}


namespace {
template <typename T, typename U>
Atop< BoundOperator< Function<T>, U >, Array >
join(
  Function<T> left,
  Partial<Operator<U>, Array> right,
  Context&
)
{
  return { { std::move(left) }, std::move(right.right) };
}
}


namespace {
template <typename T, typename U, typename V>
Partial<
  Operator<T>,
  Atop< BoundOperator< Function<U>, V >, Array >
>
join(
  Operator<T> left,
  Atop< BoundOperator< Function<U>, V >, Array > right,
  Context&
)
{
  return { std::move(left), std::move(right) };
}
}


namespace {
template <typename T, typename U, typename V>
Fork<
  Fork< Array, Operator<T>, Function<U> >,
  Operator<V>,
  Array
>
join(
  Array left,
  Partial<
    Operator<T>,
    Atop< BoundOperator< Function<U>, V >, Array >
  > right,
  Context&
)
{
  return {
    {
      std::move(left),
      std::move(right.left),
      std::move(right.right.left.left)
    },
    {},
    std::move(right.right.right)
  };
}
}


namespace {
template <typename T, typename U, typename V>
Fork<
  Fork< Values, Operator<T>, Function<U> >,
  Operator<V>,
  Array
>
join(
  Array left,
  Fork<
    Fork< Array, Operator<T>, Function<U> >,
    Operator<V>,
    Array
  > right,
  Context&
)
{
  return {
    {
      { std::move(left), std::move(right.left.left) },
      std::move(right.left.mid),
      std::move(right.left.right),
    },
    std::move(right.mid),
    std::move(right.right)
  };
}
}


namespace {
Fork<vector<Var>, Keyword<Assign>, Array>
join(
  vector<Var> left,
  Partial< Keyword<Assign>, Array > right,
  Context&
)
{
  return { std::move(left), std::move(right.left), std::move(right.right) };
}
}


namespace {
Array evaluate(Fork<Var, Keyword<Assign>, Array> arg, Context&)
{
  assert(arg.left.ptr);
  *arg.left.ptr = std::move(arg.right);
  return Array(*arg.left.ptr);
}
}


namespace {
template <typename T>
Atop<Function<T>, Array>
join(
  Function<T> left,
  Fork<Var, Keyword<Assign>, Array> right,
  Context& context
)
{
  return { std::move(left), evaluate(std::move(right), context) };
}
}


template <typename T>
static T combine(Context &, T arg)
{
  return arg;
}


template <typename Arg1, typename Arg2, typename ...Args>
static auto combine(Context &context, Arg1 arg1, Arg2 arg2, Args ...args)
{
  auto right = combine(context, std::move(arg2), std::move(args)...);
  auto left = evaluate(std::move(arg1), context);
  return join(std::move(left), std::move(right), context);
}


namespace {
Array
evaluate(Fork<vector<Var>, Keyword<Assign>, Array> arg, Context& context)
{
  if (isScalar(arg.right)) {
    assert(false);
  }

  if (arg.right.shape().size() == 1) {
    int n = arg.left.size();

    if (n != arg.right.shape()[0]) {
      assert(false);
    }

    Values values;

    for (int i=0; i!=n; ++i) {
      values.push_back(
        evaluate(
          fork(arg.left[i], arg.mid, std::move(arg.right.values()[i])),
          context
        )
      );
    }

    return Array(arg.right.shape(), std::move(values));
  }

  assert(false);
}
}


static Array posIndex(int i, const Array::Shape &shape)
{
  Values result;
  int n = shape.size();

  for (int j=n; j!=0; ) {
    --j;
    result.push_back(i % shape[j] + 1);
    i /= shape[j];
  }

  std::reverse(result.begin(), result.end());
  return makeArrayFromValues(result);
}


namespace {
Array evaluate(Atop<Function<Where>, Array> arg, Context&)
{
  if (isScalar(arg.right)) {
    assert(false);
  }

  if (arg.right.shape().size() == 1) {
    assert(false);
  }

  if (arg.right.shape().size() == 2) {
    int n = arg.right.values().size();
    Values values;

    for (int i=0; i!=n; ++i) {
      if (arg.right.values()[i] == 1) {
        values.push_back(posIndex(i, arg.right.shape()));
      }
    }

    int n2 = values.size();
    return { {n2}, std::move(values) };
  }

  assert(false);
}
}


template <typename...Args>
static auto evaluateInContext(Context &context, Args &&...args)
{
  auto x = combine(context, MakeRValue<Args>()(std::forward<Args>(args))...);
  return evaluate(std::move(x), context);
}


template <typename F>
static auto evaluateExprInContext(const F &f, Context &context)
{
  auto eval = [&](auto &&...args)
  {
    return
      evaluateInContext(
        context,
        std::forward<decltype(args)>(args)...
      );
  };

  return f(eval);
}


template <typename T>
static Array evaluateDfn(const Dfn<T> &dfn, Array arg, Context &context)
{
  Context c{context.random_engine};
  c.right_arg_ptr = &arg;
  return evaluateExprInContext(dfn.f, c);
}


namespace {
template <typename T>
Array evaluate(Atop<Dfn<T>,Array> arg, Context &context)
{
  return evaluateDfn(arg.left, std::move(arg.right), context);
}
}


namespace {
template <typename T>
auto evaluate(Dfn<T> arg, Context&)
{
  return arg;
}
}


static Array makeArray(Array a)
{
  return a;
}


static inline Array makeArray(Var a)
{
  assert(a.ptr);
  return Array(*a.ptr);
}


template <typename F>
static auto evaluateExpr(Expr<F> &&expr)
{
  auto result = evaluateExprInContext(expr.f, expr.context);
  expr.evaluated = true;
  return result;
}


namespace {
Atop<Array,Index<Array>>
join(Array left, Index<Array> right, Context &)
{
  return { std::move(left), std::move(right) };
}
}


namespace {
template <typename T>
Atop<Array,Index<Array>>
join(Var left, Index<Expr<T>> right, Context &context)
{
  Array i = makeArray(evaluateExpr(std::move(right.arg)));
  return join(Array(*left.ptr), Index<Array>{std::move(i)}, context);
}
}


namespace {
template <typename T>
Atop<Array,Index<Array>>
join(Array left, Index<Expr<T>> right, Context &context)
{
  Array i = evaluateExpr(std::move(right.arg));
  return join(std::move(left), Index<Array>{std::move(i)}, context);
}
}


template <typename F>
Expr<F>::~Expr()
{
  if (!evaluated) {
    evaluateExpr(std::move(*this));
  }

  assert(evaluated);
}


namespace {
template <typename F>
auto evaluate(Expr<F> expr, Context &context)
{
  auto result = evaluateExprInContext(expr.f, context);
  expr.evaluated = true;
  return result;
}
}


namespace {
template <typename T, typename F>
auto join(T left, Expr<F> right, Context &context)
{
  return join(std::move(left), evaluateExpr(std::move(right)), context);
}
}


namespace {
Atop<BoundOperator<Function<Plus>,Reduce>,Array>
join(
  Atop<BoundOperator<Function<Plus>, Reduce>, Function<Iota> > left,
  Array right,
  Context &context
)
{
  return {
    left.left,
    evaluate( atop(std::move(left.right), std::move(right)) , context)
  };
}
}


template <typename ...Args>
static auto defer(Args &&...args)
{
  return [&](auto f){ return f(std::forward<decltype(args)>(args)...); };
}


template <typename ...Args>
static auto makeExpr(Context &context, Args &&...args)
{
  auto f = defer(std::forward<Args>(args)...);
  return Expr<decltype(f)>{context, std::move(f)};
}


namespace {
struct Placeholder {
  std::mt19937 random_engine;
  Context context{random_engine};

  template <typename ...Args>
  auto operator()(Args &&...args)
  {
    return makeExpr(context, std::forward<decltype(args)>(args)...);
  }

  static constexpr Function<Shape>     shape = {};
  static constexpr Function<Reshape>   reshape = {};
  static constexpr Function<First>     first = {};
  static constexpr Function<Equal>     equal = {};
  static constexpr Function<Plus>      plus = {};
  static constexpr Function<Minus>     minus = {};
  static constexpr Function<Times>     times = {};
  static constexpr Function<Divide>    divide = {};
  static constexpr Function<Power>     power = {};
  static constexpr Function<Greater>   greater = {};
  static constexpr Function<Iota>      iota = {};
  static constexpr Function<Roll>      roll = {};
  static constexpr Function<Replicate> replicate = {};
  static constexpr Function<MemberOf>  member_of = {};
  static constexpr Function<Not>       isnot = {};
  static constexpr Function<Drop>      drop = {};
  static constexpr Function<Enclose>   enclose = {};
  static constexpr Function<GradeUp>   grade_up = {};
  static constexpr Function<Where>     where = {};
  static constexpr Operator<Each>      each = {};
  static constexpr Operator<Reduce>    reduce = {};
  static constexpr Operator<Product>   product = {};
  static constexpr Operator<Outer>     outer = {};
  static constexpr Operator<Beside>    beside = {};
  static constexpr Keyword<Empty>      empty = {};
  static constexpr Keyword<Assign>     assign = {};
  static constexpr Keyword<RightArg>   right_arg = {};

  template <typename ...Args>
  constexpr auto dfn(Args &&...args)
  {
    auto f = defer(std::forward<Args>(args)...);
    return Dfn<decltype(f)>{std::move(f)};
  }

  template <typename ...Args>
  constexpr auto index(Args &&...args)
  {
    auto x = makeExpr(context, std::forward<decltype(args)>(args)...);
    return Index<decltype(x)>{std::move(x)};
  }
};
}


constexpr Function<Shape>     Placeholder::shape;
constexpr Function<Reshape>   Placeholder::reshape;
constexpr Function<First>     Placeholder::first;
constexpr Function<Equal>     Placeholder::equal;
constexpr Function<Plus>      Placeholder::plus;
constexpr Function<Minus>     Placeholder::minus;
constexpr Function<Times>     Placeholder::times;
constexpr Function<Divide>    Placeholder::divide;
constexpr Function<Power>     Placeholder::power;
constexpr Function<Greater>   Placeholder::greater;
constexpr Function<Iota>      Placeholder::iota;
constexpr Function<Roll>      Placeholder::roll;
constexpr Function<Replicate> Placeholder::replicate;
constexpr Function<Drop>      Placeholder::drop;
constexpr Function<MemberOf>  Placeholder::member_of;
constexpr Function<Not>       Placeholder::isnot;
constexpr Function<Enclose>   Placeholder::enclose;
constexpr Function<GradeUp>   Placeholder::grade_up;
constexpr Function<Where>     Placeholder::where;
constexpr Keyword<Empty>      Placeholder::empty;
constexpr Keyword<Assign>     Placeholder::assign;
constexpr Keyword<RightArg>   Placeholder::right_arg;
constexpr Operator<Each>      Placeholder::each;
constexpr Operator<Reduce>    Placeholder::reduce;
constexpr Operator<Product>   Placeholder::product;
constexpr Operator<Outer>     Placeholder::outer;
constexpr Operator<Beside>    Placeholder::beside;


template <typename T>
static std::string str(const T &x)
{
  std::ostringstream stream;
  stream << x;
  return stream.str();
}


template <typename F>
static std::string str(Expr<F> expr)
{
  return str(evaluateExpr(std::move(expr)));
}


template <typename F>
Expr<F>::operator Array() &&
{
  return evaluateExpr(std::move(*this));
}


static const Array::Shape &shapeOf(const Array &a)
{
  return a.shape();
}


static void runSimpleTests()
{
  Placeholder _;
  assert(_(1) == _(1));
  assert(!(_(2) == _(3)));
  assert(_(1,2) == _(1,2));
  assert(_(1,2,3) == _(1,2,3));
  assert(shapeOf(_(1,2,3)) == Array::Shape{3});
  assert(_(_.shape, 1,2,3) == _(3));
  assert(_(3, _.equal, 3) == _(1));
  assert(_(3, _.equal, _.shape, 1,2,3) == _(1));
  assert(shapeOf(_(_.shape,1,2,3)) == Array::Shape{1});
  assert(_(1,2,3, _.equal, 1,2,3) == _(1,1,1));
  assert(_(2, _.plus, 2) == _(4));
  assert(_(4, _.minus, 1) == _(3));
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
  assert(str(_(2,2, _.reshape, 1,2,3,4)) == "( [ 1 2 ] [ 3 4 ] )");
  assert(_(1,0,2, _.replicate, 1,2,3) == _(1,3,3));
  assert(_(1, _.member_of, 1) == _(1));
  assert(_(1, _.member_of, 1,2,3) == _(1));
  assert(_(1,2,3, _.member_of, 2) == _(0,1,0));
  assert(_(_.isnot, 1) == _(0));
  assert(_(_.shape, _.shape, _.enclose, 1, 2) == _(0));
  assert(_(_.first, _.enclose, 1,2) == _(1,2));
  assert(_(5,3,9, _.index(3)) == 9);
  assert(_(_.grade_up, 5,3,9) == _(2,1,3));
  assert(_(2, _.greater, 1)== _(1));
  assert(_(2, _.power, 3) == _(8));
  assert(_(_.dfn(1, _.plus, _.right_arg), 2) == _(3));
  assert(_(_.dfn(_(_.right_arg)), 5) == _(5));
  assert(_(1, _.plus, _.beside, _.divide, 2) == 1.5);
}


static void testAssignment()
{
  Placeholder _;

  {
    Array R = 5;
    _(R, _.assign, 1, _.drop, _.iota, 5);
    assert(R == _(2,3,4,5));
  }

  {
    Array R = 3;
    assert(_(_.iota, R) == _(1,2,3));
  }

  {
    Array R = 1;
    assert(_(R, _.plus, R) == _(2));
  }

  {
    Array R = 1;
    assert(_(R, _.plus, R, _.assign, 2) == _(4));
  }

  {
    Array R = 1;
    assert(_(_(R, _.plus, R), _.plus, R, _.assign, 2) == _(6));
  }

  {
    Array R = 5;
    assert(_(R, _.assign, 1, _.drop, _.iota, R) == _(2,3,4,5));
  }
}


static void testExamples()
{
  Placeholder _;

  {
    Array result = _(_(_.iota, 2), _.outer, _.product, _.times, _.iota, 2);
    assert(result == _(2,2, _.reshape, 1, 2, 2, 4));
  }

  {
    Array R = 5;

    assert(
      _(
        _(_.isnot, R, _.member_of, R, _.outer, _.product, _.times, R),
        _.replicate, R, _.assign, 1, _.drop, _.iota, R
      ) == _(2,3,5)
    );
  }

  {
    Array X = _(1,2);

    assert(
      _(_(_.plus, _.reduce, X), _.divide, _.shape, X) == _(1, _.reshape, 1.5)
    );
  }

  {
    Array X(0);
    Array Y = _(X, _.index(_.grade_up, X, _.assign, 6, _.roll, 40));
    assert(_(_.shape, Y) == _(6));
  }

  {
    Array txt = _("<b>blah</b>");
    Array result = _(_.dfn(_.right_arg, _.member_of, "<>"), txt);
    assert(result == _(1,0,1,0,0,0,0,1,0,0,1));
  }

  {
    Array result = _(1,2,3, _.plus, _.enclose, 4,5,6);
    Array expected = _(_(5,6,7), _(6,7,8), _(7,8,9));
    assert(result == expected);
  }

  {
    Array result = _(_(5,6), _.index(2,2,_.reshape,1,2,2,1));
    assert(result == _(2,2, _.reshape, 5,6,6,5));
  }

  {
    Array W = 3;
    Array H = 2;
    Array grid = _(-1, _.plus, _.iota, H, W);
    assert(grid.shape().size() == 2);

    std::string expected =
      "( [ ( 0 0 ) ( 0 1 ) ( 0 2 ) ]"
       " [ ( 1 0 ) ( 1 1 ) ( 1 2 ) ] )";

    assert(str(grid) == expected);
  }
}


static void testCircle()
{
  Placeholder _;
  Array W = 20;
  Array H = 10;
  Array grid = _(-1, _.plus, _.iota, H, W);

  Array flags =
    _(
      _.dfn(
        .1, _.greater, _.plus, _.reduce,
        _(-.5, _.plus, _.right_arg, _.divide, H, W, _.minus, 1),
        _.power, 2
      ),
      _.each, grid
    );

  Array result = _(_(' ','*'), _.index(flags, _.plus, 1));
  std::string s;

  for (int i=0; i!=result.shape()[0]; ++i) {
    for (int j=0; j!=result.shape()[1]; ++j) {
      s += result.values()[i*result.shape()[1] + j].asCharacter();
    }

    s += '\n';
  }

  std::string expected =
    "                    \n"
    "                    \n"
    "       ******       \n"
    "     **********     \n"
    "    ************    \n"
    "    ************    \n"
    "     **********     \n"
    "       ******       \n"
    "                    \n"
    "                    \n";

  assert(s == expected);
}


static void testGrille()
{
  Placeholder _;
  Array grid = 0, grille = 0;

  _(_(grid, grille), _.assign, 5, 5, _.beside, _.reshape, _.each,
    "VRYIALCLQIFKNEVPLARKMPLFF", "XXX X XXX X X XXX XXX  XX");

  Array result = _(grid, _.index(_(_.where, grille, _.equal, ' ')));
  assert(result == makeArrayFromString("ILIKEAPL"));
}


int main()
{
  runSimpleTests();
  testAssignment();
  testExamples();
  testCircle();
  testGrille();
}
