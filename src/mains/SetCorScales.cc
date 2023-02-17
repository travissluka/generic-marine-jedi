
/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "genericMarine/Traits.h"

#include "oops/runs/Run.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace genericMarine {

class SetCorScales : public oops::Application {
 public:
  explicit SetCorScales(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}
  static const std::string classname() { return "genericMarine::SetCorScales";}

  int execute(const eckit::Configuration & config, bool /*validate*/) const {
    // setup geometry
    const eckit::LocalConfiguration geometryConfig(config, "geometry");
    const Geometry geom(geometryConfig, this->getComm());

    // setup increment
    const util::DateTime date;
    const oops::Variables vars(config, "variables");
    ASSERT(vars.size() == 1);  // TODO(travis)  make work with >1 field
    Increment rh(geom, vars, date);

    // get rh parameters
    // rh is calculated as follows :
    // 1) rh = "base value" + rossby_radius * "rossby mult"
    // 2) minimum value of "min grid mult" * grid_size is imposed
    // 3) min/max are imposed based on "min value" and "max value"
    // 4) converted from a gaussian sigma to Gaspari-Cohn cutoff distance
    const eckit::LocalConfiguration corrConf(config, "scales");
    double baseValue = corrConf.getDouble("base value", 0.0);
    double rossbyMult = corrConf.getDouble("rossby mult", 0.0);
    double minGridMult = corrConf.getDouble("min grid mult", 0.0);
    double minValue = corrConf.getDouble("min value", 0.0);
    double maxValue = corrConf.getDouble("max value",
                                       std::numeric_limits<double>::max());

    // get required fields from geometry
    auto rossbyRadius = atlas::array::make_view<double, 2>(
      geom.extraFields().field("rossby_radius"));
    auto area = atlas::array::make_view<double, 2>(
        geom.extraFields().field("area"));

    // create the atlas field
    atlas::FieldSet param_fieldSet;
    atlas::Field param_field = geom.functionSpace().createField<double>(
      atlas::option::levels(1) | atlas::option::name(vars[0]));
    param_fieldSet.add(param_field);
    auto param_view = atlas::array::make_view<double, 2>(param_field);

    // start assigning values to the field
    param_view.assign(baseValue);

    for ( int i = 0; i < param_field.size(); i++ ) {
      param_view(i, 0) += rossbyMult * rossbyRadius(i, 0);
      param_view(i, 0) = std::max(param_view(i, 0),
                                  minGridMult*sqrt(area(i, 0)));
      param_view(i, 0) = std::max(param_view(i, 0), minValue);
      param_view(i, 0) = std::min(param_view(i, 0), maxValue);

      // note: BUMP expects the length as a Gaspari-Cohn cutoff length,
      //   but we probably think of it as a Gaussian 1 sigma, so convert.
      param_view(i, 0) *= 3.57;  // gaussian to GC factor
    }

    // save rh
    rh.fromFieldSet(param_fieldSet);
    oops::Log::test() << "Output horizontal scales: " << rh << std::endl;
    const eckit::LocalConfiguration rhOutputConfig(config, "output.rh");
    rh.write(rhOutputConfig);

    // generate/save rv
    const eckit::LocalConfiguration rvOutputConfig(config, "output.rv");
    Increment rv(rh);
    param_view.assign(1.0);
    rv.fromFieldSet(param_fieldSet);
    oops::Log::test() << "Output vertical scales: " << rv << std::endl;
    rv.write(rvOutputConfig);

    return 0;
  }

 private:
  std::string appname() const { return "genericMarine::SetCorScales";}
};

}  // namespace genericMarine
// -----------------------------------------------------------------------------

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  genericMarine::SetCorScales setCorScales;
  return run.execute(setCorScales);
}
