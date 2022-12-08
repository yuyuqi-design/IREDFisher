from enzyme_evolver.app.app import create_app
import os

production_mode = os.environ.get('PRODUCTION') or False

print('Creating app')
main_app = create_app(use_talisman=production_mode)

if __name__ == '__main__':
    main_app.run(debug=not production_mode)
