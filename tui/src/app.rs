use crate::modules::{
    math_eval::MathEvalModule, menu::MenuModule, newton::NewtonModule,
    placeholder::PlaceholderModule, Modules, MODULES,
};

use std::time::Duration;

use color_eyre::Result;
use crossterm::event::KeyEvent;
use ratatui::{
    crossterm::event::{Event, EventStream, KeyCode, KeyEventKind},
    layout::{Constraint, Layout},
    style::Stylize,
    text::Line,
    DefaultTerminal, Frame,
};
use tokio_stream::StreamExt;

#[derive(Debug, Default, PartialEq, Eq)]
enum State {
    #[default]
    Running,
    Quit,
}

#[derive(Debug, Default)]
enum Mode {
    #[default]
    Normal,
    Insert,
}

#[derive(Debug, Default)]
pub struct App {
    current_state: State,
    current_mode: Mode,
    current_module: Modules,
    menu_module: MenuModule,
    math_eval_module: MathEvalModule,
    newton_module: NewtonModule,
    placeholder_module: PlaceholderModule,
}

impl App {
    const FRAMES_PER_SECOND: f32 = 60.0;
    const MAPPINGS: &str = "q: Quit | i/a: Insert | j/k/Down/Up: Move | Enter: Run | Tab: Next Module | BackTab: Prev Module | Esc: End Input";

    pub async fn run(mut self, mut terminal: DefaultTerminal) -> Result<()> {
        let period = Duration::from_secs_f32(1.0 / Self::FRAMES_PER_SECOND);
        let mut interval = tokio::time::interval(period);
        let mut events = EventStream::new();

        // Main event loop
        while self.current_state != State::Quit {
            tokio::select! {
                _ = interval.tick() => { terminal.draw(|frame| self.draw(frame))?; },
                Some(Ok(event)) = events.next() => self.handle_event(&event),
            }
        }
        Ok(())
    }

    fn draw(&self, frame: &mut Frame) {
        // Setup the general layout
        let vertical = Layout::vertical([
            Constraint::Length(1),
            Constraint::Fill(1),
            Constraint::Length(1),
        ]);
        let [title_area, body_area, info_area] = vertical.areas(frame.area());

        // Setup title widget
        let title = Line::from(self.current_module.to_string())
            .centered()
            .bold();

        // Setup info layout
        let info_horizontal =
            Layout::horizontal([Constraint::Percentage(95), Constraint::Length(20)]);
        let [info_mappings_area, info_mode_area] = info_horizontal.areas(info_area);

        // Setup info widgets
        let info_mappings = Line::from(Self::MAPPINGS).left_aligned();
        let info_mode = match self.current_mode {
            Mode::Normal => Line::from("Normal Mode").blue(),
            Mode::Insert => Line::from("Insert Mode").green(),
        }
        .right_aligned();

        // Setup body layout
        let horizontal =
            Layout::horizontal([Constraint::Percentage(20), Constraint::Percentage(80)]);
        let [body_area_left, body_area_right] = horizontal.areas(body_area);

        // Render everything
        frame.render_widget(title, title_area);
        frame.render_widget(&self.menu_module, body_area_left);
        match self.current_module {
            Modules::ExprEval => frame.render_widget(&self.math_eval_module, body_area_right),
            Modules::Newton => frame.render_widget(&self.newton_module, body_area_right),
            Modules::Placeholder => frame.render_widget(&self.placeholder_module, body_area_right),
        }
        frame.render_widget(info_mappings, info_mappings_area);
        frame.render_widget(info_mode, info_mode_area);
    }

    fn handle_event(&mut self, event: &Event) {
        if let Event::Key(key) = event {
            if key.kind == KeyEventKind::Press {
                match self.current_mode {
                    Mode::Normal => match key.code {
                        KeyCode::Char('q') => self.current_state = State::Quit,
                        KeyCode::Char('i') | KeyCode::Char('a') => self.begin_input(),
                        KeyCode::Char('j') | KeyCode::Down => self.prev_field(),
                        KeyCode::Char('k') | KeyCode::Up => self.next_field(),
                        KeyCode::Enter => self.module_run(),
                        KeyCode::Tab => {
                            self.menu_module.next_module();
                            let index = self.menu_module.get_selected_module();
                            self.current_module = MODULES[index];
                        }
                        KeyCode::BackTab => {
                            self.menu_module.prev_module();
                            let index = self.menu_module.get_selected_module();
                            self.current_module = MODULES[index];
                        }
                        _ => {}
                    },
                    Mode::Insert => match key.code {
                        KeyCode::Tab => {
                            self.end_input();
                            self.menu_module.next_module();
                            let index = self.menu_module.get_selected_module();
                            self.current_module = MODULES[index];
                        }
                        KeyCode::BackTab => {
                            self.end_input();
                            self.menu_module.prev_module();
                            let index = self.menu_module.get_selected_module();
                            self.current_module = MODULES[index];
                        }
                        KeyCode::Esc | KeyCode::Enter => self.end_input(),
                        _ => self.input(*key),
                    },
                }
            }
        }
    }

    fn begin_input(&mut self) {
        self.current_mode = Mode::Insert;
        match self.current_module {
            Modules::ExprEval => self.math_eval_module.insert_mode_enter(),
            Modules::Newton => self.newton_module.insert_mode_enter(),
            Modules::Placeholder => {}
        }
    }

    fn end_input(&mut self) {
        self.current_mode = Mode::Normal;
        match self.current_module {
            Modules::ExprEval => self.math_eval_module.insert_mode_exit(),
            Modules::Newton => self.newton_module.insert_mode_exit(),
            Modules::Placeholder => {}
        }
    }

    fn input(&mut self, c: KeyEvent) {
        match self.current_module {
            Modules::ExprEval => self.math_eval_module.input(c),
            Modules::Newton => self.newton_module.input(c),
            Modules::Placeholder => {}
        }
    }

    fn prev_field(&mut self) {
        match self.current_module {
            Modules::ExprEval => self.math_eval_module.prev_field(),
            Modules::Newton => self.newton_module.prev_field(),
            Modules::Placeholder => {}
        }
    }

    fn next_field(&mut self) {
        match self.current_module {
            Modules::ExprEval => self.math_eval_module.next_field(),
            Modules::Newton => self.newton_module.next_field(),
            Modules::Placeholder => {}
        }
    }

    fn module_run(&mut self) {
        match self.current_module {
            Modules::ExprEval => self.math_eval_module.run(),
            Modules::Newton => self.newton_module.run(),
            Modules::Placeholder => {}
        }
    }
}
