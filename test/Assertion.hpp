class Assertion {
   private:
    int err;
    int last_err;

   public:
    Assertion() : err(0) {}

    void operator()(bool t) {
        err = err == 0 && t ? (last_err = 0) : (last_err = 1);
    }
    int last_status() const { return last_err; }
    int status() const { return err; }
};