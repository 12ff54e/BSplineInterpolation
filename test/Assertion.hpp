class Assertion {
   private:
    int err;

   public:
    Assertion() : err(0) {}

    void operator()(bool t) { err = err == 0 && t ? 0 : 1; }
    int status() const { return err; }
};